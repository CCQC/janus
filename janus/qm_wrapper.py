"""
QM wrapper super class
"""
class QM_wrapper(object):

    def __init__(self, system):

        self._system = system
        self._qm = None
    
    def qm_info(self):
        pass

    def get_qm(self):
        if self._qm is None:
            self.qm_info()
        return self._qm 
