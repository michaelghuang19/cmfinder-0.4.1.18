
class MarkovModelStats {
	int order;
	int alphabetSize;
	VariableDimVector<double> data;
	std::list<int> path;

	MarkovModelStats (void);

	template <class Container>
	double GetProbOfNuc_templ (int nucForPr,const Container& container);
public:
	MarkovModelStats (FILE *file); // calls LoadFromDump
	MarkovModelStats (int order_,int alphabetSize_=4);
	~MarkovModelStats ();
	void ClearData (void);
	void PushOnPath (int nuc);
	void Dump (FILE *out);
	void LoadFromDump (FILE *file);
	void LoadFromDump (const char *fileName);
	// get raw number
	double Get (void);
	void Set (double t);
	void Inc (void); // equiv to Set(Get()+1)
	int GetOrder (void);
	// get Pr(nuc|context)
	double GetProbOfNuc (int nucForPr,const std::list<int>& context); // context is a list of 'order' nucs, with context.back() being the most recent
	double GetProbOfNuc (int nucForPr,const NaryCounter& context);
	double GetProbOfNuc (int nucForPr,const vector<int>& context);
	// convenience function for 0-order models
	double GetProbOfNuc_0order (int nucForPr);

	// a Markov model assigns Pr(x|Y) where x is a nuc, and Y is a string of nucs of length 'order'.  I'm calling Y the 'context'.  In some cases, we'd like to find the probability of these contexts according to the model, e.g. to initialize the context for generating a random string according to the model.
	void GetContextDistribution (VariableDimVector<double>& probabilityOfEachContext);

	// makes 0-order Markov model with uniform stats
	static MarkovModelStats *NewUniformMarkov (void);
	// makes 0-order Markov model with given nucleotide probs (nucProbs is array of size 4)
	static MarkovModelStats *NewMarkov0 (double *nucProbs);
	// makes a Markov model of order (order-1), where 'order' is the order of the given MM.  This is pretty much just for testing my code
	static MarkovModelStats *NewDecrementedOrderMarkov (MarkovModelStats *inputMM);
};
