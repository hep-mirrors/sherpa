from os import path, makedirs

class test_process(object):

    def __init__(self, in_ids, out_ids, model, threshold=1.0e-9, cms=2000.0):
        self._in_ids        = in_ids
        self._out_ids       = out_ids
        self._model         = model
        self._threshold     = threshold
        self._cms           = cms

    def test_dir_name(self):
        return '_'.join([self._model._name]+[str(id) for id in self._in_ids]+[str(id) for id in self._out_ids])

    def __str__(self):
        return ' '.join([str(ide) for ide in self._in_ids])+' -> '+' '.join([str(ide) for ide in self._out_ids])+' ({0})'.format(self._model._name)
