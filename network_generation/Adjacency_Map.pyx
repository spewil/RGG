#distutils: language = c++
#cython: boundscheck = False
STUFF = "HI"

# Import the map and vector templates from the STL
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from libcpp.utility cimport pair as cpp_pair

ctypedef cpp_vector[int] cpp_neighbourhood
ctypedef cpp_map[int, cpp_neighbourhood] cpp_adjacency_map
ctypedef cpp_pair[int, cpp_neighbourhood] cpp_item

# Import a few operators because they aren't supported by cython syntax
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc


#Define adjacency map class for detecting the giant component:
cdef class AdjacencyMap:
    
    cdef:
        cpp_adjacency_map container
        
    def __init__(self, int[:, :] adjacency_list):
        cdef:
            int i, ego, alter
            cpp_neighbourhood neighbourhood
            
        # Iterate over all entries of the adjacency list
        for i in range(adjacency_list.shape[0]):
            ego = adjacency_list[i, 0]
            alter = adjacency_list[i, 1]
            
            # Check if the ego is already in the map 
            # (see http://stackoverflow.com/a/101980/1150961 for details)
            lb = self.container.lower_bound(ego)
            
            # Check if the key already exists
            if lb != self.container.end() and ego == deref(lb).first:
                # Add the node to the pair
                deref(lb).second.push_back(alter)
            else:
                # Insert a new key value pair
                neighbourhood = cpp_neighbourhood()
                neighbourhood.push_back(alter)
                self.container.insert(lb, cpp_item(ego, neighbourhood))
                
    def get(self, int ego):
        """
        Get the neighbours of `ego` or `None` if `ego` isn't in the map.
        """
        # Search the dictionary
        iterator = self.container.find(ego)
        # Return none if we didn't find anything
        if iterator == self.container.end():
            return None
        # Create a list of values from the vector
        values = []
        # Iterate over the neighbourhood and add values to the list
        neighbourhood = deref(iterator).second
        neighbourhood_iterator = neighbourhood.begin()
        while neighbourhood_iterator != neighbourhood.end():
            values.append(deref(neighbourhood_iterator))
            preinc(neighbourhood_iterator)
            
        return values
    
    def _get_many(self, int ego, int repeats):
        """
        Simple function to illustrate overhead.
        """
        cdef int i
        # Try to find the ego a large number of times
        for i in range(repeats):
            iterator = self.container.find(ego)


