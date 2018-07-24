"""
AQMMM class for adaptive QMMM computations
"""
class AQMMM(object):

    def __init__(self, scheme, partition_scheme):

        self.scheme = scheme
        self.partition_scheme = partition_scheme

    def partition(self):

        self.define_buffer_zone()

        # make this into class structure?
        if self.scheme == 'ONIOM-XS':
            self.oniom_xs(partition=True)

    def save(self):
        pass

    def get_info(self):

        if self.scheme == 'ONIOM-XS':
            self.oniom_xs(get_info=True)
        
        return self.info

    def define_buffer_zone(self):

        if self.partition_scheme == 'distance': 
            pass
            

    def oniom_xs(partition=False, get_info=False):
        
        if partition is True:
            pass

        if get_info is True:
            pass


