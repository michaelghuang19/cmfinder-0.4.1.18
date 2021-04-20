/*
NaryCounter:
Counts in base-N with D digits, i.e. exploring all combinations.
Assumes that N<256, since otherwise going thru all combinations would take too long

  See also 'BooleanCounter.h'

  example code:

  NaryCounter counter(numDigits,base);
  counter.Init();
  bool counting=true;
  while (counting) {

	// use counter[0...(numDigits-1)]

	counting=counter.Next();
  }
*/

class NaryCounter {
protected:
	vector<char> array;
	int numDigits;
	char base;
public:
	NaryCounter (int _numDigits,char _base);
	~NaryCounter ();

	// re-initialize to all 0s
	void Init (void);

	// increment by 1
	bool /* has next */ Next (void);

	inline int operator [] (int i) const {
		return array[i];
	}

	// for convenience, make this look like a container
	typedef vector<char>::const_iterator const_iterator;
	const_iterator begin (void) const;
	const_iterator end (void) const;
};
