
cdef extern from "Model/Sequence/Data.h" namespace "GPPG::Model":
	cdef cppclass SequenceData:
		SequenceData(int) except +
		int length()
		short get(int i)
		void set(int i, short c)
		
cdef extern from "Model/Sequence/Operation.h" namespace "GPPG::Model":
	cdef cppclass SequenceRoot:
		SequenceRoot( SequenceData* ) except +

				
cdef class Sequence:
	cdef SequenceData *thisptr
		
	cpdef int length(self)
		
	cpdef short get(self, int i)
		
	cpdef set(self, int i, short c)

	