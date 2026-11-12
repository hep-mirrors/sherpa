"""Watches the Sherpa jobs of run_mb.py and stops the ones that are stuck or too slow, so that a
round ends by a deadline.  A stopped job keeps the events it has made: Sherpa is stopped, the
reader (mb_analysis.py) sees the end of the pipe and writes out.json as usual, and the file
stopped_early in the job folder says why.  The fit then uses fewer events for that point, with
correspondingly larger MC errors.

    python3 watch_jobs.py <run_mb.py round> <deadline HH:MM>     e.g.  round5 23:50

Every minute, for every running Sherpa:
  stuck      run.log (Sherpa's progress line) has not changed for STUCK_MINUTES      -> stop
  too slow   running longer than SLOW_MINUTES and at least MIN_FRACTION of its events -> stop
At the deadline: run_mb.py is stopped first (no new jobs), then every running Sherpa.
The watcher ends when no Sherpa of the round is running and run_mb.py has ended.
"""
import os
import re
import signal
import subprocess
import sys
import time

from paths import TUNING

SHERPA = "/Users/ekremorhandemiray/code_software/sherpaMainBranch/build_lund/bin/Sherpa"
STUCK_MINUTES = 30
SLOW_MINUTES = 60
MIN_FRACTION = 0.2


def log(text):
    print(time.strftime("%H:%M"), text, flush=True)


def running_sherpas():
    "[(pid, folder, seconds running)] of the Sherpa processes (not their caffeinate parents)"
    output = subprocess.run(["ps", "-axo", "pid=,etime=,command="], capture_output=True, text=True).stdout
    result = []
    for line in output.split("\n"):
        columns = line.split(None, 2)
        if len(columns) < 3:
            continue
        if not columns[2].startswith(SHERPA):
            continue
        pid = int(columns[0])
        lsof = subprocess.run(["lsof", "-a", "-p", str(pid), "-d", "cwd", "-Fn"], capture_output=True, text=True).stdout
        folder = None
        for entry in lsof.split("\n"):
            if entry.startswith("n"):
                folder = entry[1:]
        if folder is None:
            continue
        result.append((pid, folder, elapsed_seconds(columns[1])))
    return result


def elapsed_seconds(etime):
    "ps etime [[dd-]hh:]mm:ss -> seconds"
    days = 0
    if "-" in etime:
        (day_part, etime) = etime.split("-")
        days = int(day_part)
    parts = etime.split(":")
    seconds = 0
    for part in parts:
        seconds = seconds * 60 + int(part)
    return days * 86400 + seconds


def events_made(folder):
    "the last 'Event N (' in run.log"
    text = open(os.path.join(folder, "run.log"), errors="replace").read()
    numbers = re.findall(r"Event (\d+) \(", text)
    if len(numbers) == 0:
        return 0
    return int(numbers[-1])


def events_wanted(folder):
    "EVENTS of card.yaml, e.g. 100k -> 100000"
    for line in open(os.path.join(folder, "card.yaml")):
        if line.startswith("EVENTS:"):
            value = line.split(":")[1].strip().strip("'\"")
            if value.endswith("k"):
                return int(float(value[:-1]) * 1000)
            return int(float(value))
    return 0


def stop(pid, folder, reason, round_name):
    """Stop one Sherpa and keep its events.  run_mb.py kills the reader of a Sherpa that ended
    with an error code (run_mb.py, run()), maybe before the reader has written out.json, so
    run_mb.py is paused (SIGSTOP) until out.json is there, at most 2 minutes."""
    made = events_made(folder)
    with open(os.path.join(folder, "stopped_early"), "w") as f:
        f.write(reason + "; " + str(made) + " of " + str(events_wanted(folder)) + " events\n")
    paused = run_mb_pids(round_name)
    for p in paused:
        os.kill(p, signal.SIGSTOP)
    os.kill(pid, signal.SIGTERM)
    waited = 0
    while not os.path.exists(os.path.join(folder, "out.json")) and waited < 120:
        time.sleep(2)
        waited += 2
    for p in paused:
        os.kill(p, signal.SIGCONT)
    if os.path.exists(os.path.join(folder, "out.json")):
        kept = "out.json written"
    else:
        kept = "NO out.json"
    log("stopped " + os.path.relpath(folder, TUNING) + ": " + reason + ", " + str(made) + " events, " + kept)


def run_mb_pids(round_name):
    output = subprocess.run(["pgrep", "-f", "run_mb.py run " + round_name], capture_output=True, text=True).stdout
    pids = []
    for line in output.split():
        pids.append(int(line))
    return pids


def main(round_name, deadline):
    round_folder = os.path.join(TUNING, "runs", round_name)
    # the deadline as a time: today at HH:MM, or tomorrow if that has passed
    (hour, minute) = deadline.split(":")
    now = time.localtime()
    deadline_time = time.mktime((now.tm_year, now.tm_mon, now.tm_mday, int(hour), int(minute), 0, 0, 0, -1))
    if deadline_time < time.time():
        deadline_time += 86400
    log("watching " + round_name + ", deadline " + time.strftime("%a %H:%M", time.localtime(deadline_time)))
    past_deadline = False
    while True:
        if not past_deadline and time.time() >= deadline_time:
            past_deadline = True
            for pid in run_mb_pids(round_name):
                os.kill(pid, signal.SIGTERM)
                log("deadline: stopped run_mb.py (pid " + str(pid) + ")")
        jobs = []
        for (pid, folder, seconds) in running_sherpas():
            if folder.startswith(round_folder):
                jobs.append((pid, folder, seconds))
        for (pid, folder, seconds) in jobs:
            made = events_made(folder)
            wanted = events_wanted(folder)
            quiet = time.time() - os.path.getmtime(os.path.join(folder, "run.log"))
            if past_deadline:
                stop(pid, folder, "deadline " + deadline, round_name)
            elif made > 0 and quiet > STUCK_MINUTES * 60:
                stop(pid, folder, "stuck: no progress for %d min" % (quiet / 60), round_name)
            elif seconds > SLOW_MINUTES * 60 and wanted > 0 and made >= MIN_FRACTION * wanted:
                stop(pid, folder, "too slow: running %d min" % (seconds / 60), round_name)
        if len(jobs) == 0 and len(run_mb_pids(round_name)) == 0:
            log("no jobs left, watcher ends")
            return
        time.sleep(60)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
