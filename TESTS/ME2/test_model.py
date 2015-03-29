from os import path
import imp

class test_model(object):

    def __init__(self, ufo_name, ufo_path, param_card_path):
        self._ufo_path   = ufo_path
        self._name       = ufo_name
        self._param_card = param_card_path
        try:
            f,pathname, desc = imp.find_module(self._name, [ufo_path])
            self._ufo_model  = imp.load_module(self._name, f, pathname, desc)
        except ImportError:
            raise ValueError("Could not import model {0} from directory {1}".format(ufo_name, ufo_path))
        
