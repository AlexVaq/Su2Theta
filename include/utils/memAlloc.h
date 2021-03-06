#ifndef	MEM_ALLOC
	#define	MEM_ALLOC

	#include "enumFields.h"

	void	trackFree    (void *);
	void	trackAlloc   (void **, size_t);
	void	alignAlloc   (void **, size_t, size_t);
	void	printMemStats();
#endif

