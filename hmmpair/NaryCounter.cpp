#include "hmmpair.h"

NaryCounter::NaryCounter (int _numDigits,char _base)
{
	numDigits=_numDigits;
	base=_base;

	Init();
}
NaryCounter::~NaryCounter ()
{
}
void NaryCounter::Init (void)
{
	array.assign(numDigits,0);
}
bool NaryCounter::Next (void)
{
	int i;
	i=0;
	while (i<numDigits) {
		array[i]++;
		if (array[i]==base) {
			array[i]=0;
			i++;
		}
		else {
			return true;
		}
	}

	return false;
}
NaryCounter::const_iterator NaryCounter::begin (void) const
{
	return array.begin();
}
NaryCounter::const_iterator NaryCounter::end (void) const
{
	return array.end();
}
