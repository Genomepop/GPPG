"""
cdef extern from "Model/Sequence/Data.h" namespace "GPPG::Model":
	cdef cppclass SequenceData:
		SequenceData(int) except +
		int length()
		short get(int i)
		void set(int i, short c)
"""		
cdef class Sequence:
	#cdef SequenceData* thisptr
	
	def __cinit__(self, int l):
		self.thisptr = new SequenceData(l)
		
	def __dealloc__(self):
		del self.thisptr
		
	cpdef int length(self):
		return self.thisptr.length()
		
	cpdef short get(self, int i):
		return self.thisptr.get(i)
		
	cpdef set(self, int i, short c):
		self.thisptr.set(i, c)
		
cdef class PySequenceRoot:
	cdef SequenceRoot* thisptr
	
	def __cinit__(self, Sequence data):
		self.thisptr = new SequenceRoot( data.thisptr )
	
