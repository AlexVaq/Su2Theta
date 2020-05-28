#ifndef	__QCHARGECLASS
	#define	__QCHARGECLASS

	#include <enumFields.h>
	#include "utils/utils.h"
	#include "simd/simd.h"

	using namespace Su2Enum;

	template<class T>
	int?	simplex (Su2 &y1, Su2 &y2, Su2 &y3, Su2 &y4) {

		T::data	d1 = y2.a[1]*y1.a[2] - y1.a[1]*y2.a[2];
		T::data	d2 = y3.a[1]*y4.a[2] - y4.a[1]*y3.a[2];
		T::data	d3 = y4.a[1]*y2.a[2] - y2.a[1]*y4.a[2];
		T::data	d4 = y2.a[1]*y3.a[2] - y3.a[1]*y2.a[2];
		T::data	d5 = y3.a[1]*y1.a[2] - y1.a[1]*y3.a[2];
		T::data	d6 = y4.a[1]*y1.a[2] - y1.a[1]*y4.a[2];

		T::data	t1 = y2.a[0]*d2 + y3.a[0]*d3 + y4.a[0]*d4;
		T::data	t2 = y4.a[0]*d5 - y3.a[0]*d6 - y1.a[0]*d2;
		T::data	t3 = y2.a[0]*d6 - y4.a[0]*d1 - y1.a[0]*d3;
		T::data	t4 = y3.a[0]*d1 - y2.a[0]*d5 - y1.a[0]*d4;

		T::data	t0 = y1.a[3]*t1 + y2.a[3]*t2 + y3.a[3]*t3 + y4.a[3]*t4;

		T::data	i1 = t0*t1;
		T::data	i2 = t0*t2;
		T::data	i3 = t0*t3;
		T::data	i4 = t0*t4;

		/* TODO

		   Check if the i's are positive. If so, add
		   orientation (MISSING!) x IPS (MISSING!) x Sign t0
		   to the topological charge
		*/
	}

