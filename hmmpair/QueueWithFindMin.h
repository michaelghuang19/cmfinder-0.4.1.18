
template <class T,class Predicate>
class PriorityHeapWithPointers {
protected:
	typedef int HeapSlot;
	Predicate predicate;
public:
	class SingletonPointer {
		friend class PriorityHeapWithPointers;
		HeapSlot heapSlot;
		void SetHeapSlot (HeapSlot heapSlot_) {
			heapSlot=heapSlot_;
		}
	public:
		SingletonPointer () {
			heapSlot=-1;
		}
		~SingletonPointer () {
		}
	};
protected:
	struct Datum {
		T key;
		SingletonPointer *pointer;
	};
	typedef vector<Datum> HeapArray;
	HeapArray heapArray;
	int maxSize,size;

	inline HeapSlot GetRoot (void) {
		return 0;
	}
	inline HeapSlot GetParent (HeapSlot heapSlot) {
		assert(heapSlot!=GetRoot());
		return heapSlot>>1;
	}
	inline HeapSlot GetLeftChild (HeapSlot heapSlot) {
		return heapSlot<<1;
	}
	inline HeapSlot GetRightChild (HeapSlot heapSlot) {
		return (heapSlot<<1)+1;
	}

	void BruteForceSwap (HeapSlot slot1,HeapSlot slot2) { // we could probably save an extra SetHeapSort call by only doing it when the element is at its final resting place... but, this'd be more difficult
		std::swap(heapArray[slot1],heapArray[slot2]);
		heapArray[slot1].pointer->SetHeapSlot(slot1); 
		heapArray[slot2].pointer->SetHeapSlot(slot2);
	}
	void PercolateUp (HeapSlot thisSlot) {
		while (1) {
			if (thisSlot==GetRoot()) {
				// we're done
				break;
			}
			else {
				HeapSlot parentSlot=GetParent(thisSlot);
				if (!(predicate(heapArray[thisSlot].key,heapArray[parentSlot].key))) {
					// alles unter Controlle
					break;
				}
				else {
					// need to percolate up
					BruteForceSwap(parentSlot,thisSlot);
					thisSlot=parentSlot; // loop again
				}
				thisSlot=parentSlot;
			}
		}
	}
	void PercolateDown (HeapSlot thisSlot) {
		while (1) {
			HeapSlot leftChild=GetLeftChild(thisSlot);
			if (leftChild>=size) {
				// we're done
				break;
			}
			HeapSlot rightChild=GetRightChild(thisSlot);
			assert(rightChild<=size); // rightChild==leftChild+1
			if (rightChild==size) { // edge case, only has left child, AND the left child is a leaf
				if (predicate(heapArray[leftChild].key,heapArray[thisSlot].key)) {
					BruteForceSwap(leftChild,thisSlot);
				}
				break; // left child is a leaf
			}
			else {
				// has two children
				if (predicate(heapArray[leftChild].key,heapArray[thisSlot].key)) {
					if (predicate(heapArray[rightChild].key,heapArray[leftChild].key)) {
						BruteForceSwap(rightChild,thisSlot);
						thisSlot=rightChild;
					}
					else {
						BruteForceSwap(leftChild,thisSlot);
						thisSlot=leftChild;
					}
				}
				else {
					if (predicate(heapArray[rightChild].key,heapArray[thisSlot].key)) {
						BruteForceSwap(rightChild,thisSlot);
						thisSlot=rightChild;
					}
					else {
						// everything's in order
						break;
					}
				}
			}
		}
	}
public:
	PriorityHeapWithPointers (int maxSize_) {
		maxSize=maxSize_;
		heapArray.resize(maxSize);
		size=0;
	}
	~PriorityHeapWithPointers () {
	}
	void reserve (int maxSize_) {
		maxSize=maxSize_;
		heapArray.resize(maxSize_);
	}

	// only 1 instance of a singleton pointer to the same referent may be live at a time.  (The issue is that pointers have to be updated when things are moved; if we have to notify more than one pointer & pointers can be deleted at any time, we'd need a doubly-linked list, which is kind of expensive
	void Enter (SingletonPointer& pointer,const T& key) {
		assert(size<maxSize);
		heapArray[size].key=key;
		heapArray[size].pointer=&pointer;
		pointer.SetHeapSlot(size);
		size++;
		PercolateUp(size-1);
	}
	void Delete (const SingletonPointer& pointer) {
		assert(heapArray[pointer.heapSlot].pointer==&pointer); // two-way pointers should be consistent

		// alg: take last element of heap & slap it into the empty slot, then percolate up or down
		if (pointer.heapSlot==size-1) {
			// easy to delete
			size--;
		}
		else {
			// move last element of heap into empty slot
			heapArray[pointer.heapSlot]=heapArray[size-1];
			heapArray[pointer.heapSlot].pointer->SetHeapSlot(pointer.heapSlot);
			size--;

			// since empty slot was greater than its parent and less than its children, the new value in empty slot can be either less than its parent, or greater than a child, but not both.  So if we call PercolateUp, then PercolateDown, we should be done (and at most one of these functions will modify something)
			PercolateUp(pointer.heapSlot);
			PercolateDown(pointer.heapSlot);
		}
	}
	const T& FindMin (void) const {
		assert(size>0);
		return heapArray[0].key;
	}
	void ClearAll (void) {
		size=0;
	}

	void VerifySanity (void) {
		for (HeapSlot slot=GetRoot(); slot!=size; slot++) {

			// verify pointer points at me
			assert(heapArray[slot].pointer->heapSlot==slot);

			// verify heap's partial order property wrt children
			HeapSlot leftChild=GetLeftChild(slot);
			if (leftChild<size) {
				assert(predicate(heapArray[slot].key,heapArray[leftChild].key));
			}
			HeapSlot rightChild=GetRightChild(slot);
			if (rightChild<size) {
				assert(predicate(heapArray[slot].key,heapArray[rightChild].key));
			}
		}
	}
};

template <class T>
class CircularArray {
protected:
	int circularFirst,circularLast;
	int size;
	vector<T> circularArray;
public:
	CircularArray (int maxSize)
	{
		size=maxSize+1;
		circularArray.resize(size);
		circularFirst=circularLast=0;
	}
	~CircularArray () {
	}
	void resize (int maxSize) {
		size=maxSize+1;
		circularArray.resize(size);
	}
	T& EnterNull () {
		circularFirst--;
		if (circularFirst<0) {
			circularFirst=size-1;
		}
		assert(circularFirst!=circularLast); // or ran out of space

		return circularArray[circularFirst];
	}
	void Enter (const T& key) {
		T& slot=EnterNull();
		slot=key;
	}
	T& LeaveNull (void) {
		assert(circularFirst!=circularLast); // or nothing to delete
		circularLast--;
		if (circularLast<0) {
			circularLast=size-1;
		}
		return circularArray[circularLast];
	}
	T& Leave (void) {
		LeaveNull();
	}
	void ClearAll (void) {
		circularFirst=circularLast=0;
	}
};

template <class T,class Predicate=std::less<T> >
class QueueWithFindMin_CircularArrayAndHeap {
protected:
	struct Datum {
		typename PriorityHeapWithPointers<T,Predicate>::SingletonPointer pointer;
	};
	PriorityHeapWithPointers<T,Predicate> priorityHeap;
	CircularArray<Datum> circularArray;
public:
	QueueWithFindMin_CircularArrayAndHeap (int maxSize) 
		: priorityHeap(maxSize),circularArray(maxSize)
	{
	}
	~QueueWithFindMin_CircularArrayAndHeap () {
	}
	void resize (int maxSize) {
		priorityHeap.reserve(maxSize);
		circularArray.resize(maxSize);
	}

	void Enter (const T& key) {
		Datum& slot=circularArray.EnterNull();
		priorityHeap.Enter(slot.pointer,key);
	}
	void Leave (void) {
		Datum& slot=circularArray.LeaveNull();
		priorityHeap.Delete(slot.pointer);
	}
	const T& FindMin (void) const {
		return priorityHeap.FindMin();
	}
	void ClearAll (void) {
		priorityHeap.ClearAll();
		circularArray.ClearAll();
	}
};

// when you enter something into a queue, and the only read operation you can perform is 'FindMin', then any item with lower priority than the entered item is effectively invisible, and can be removed.  I'd like to explore the practical impact of this idea, so I'm using the STL multiset, which is an easy implementation
template <class T,class Predicate=std::less<T> >
class OpaqueQueueWithFindMin_CircularArrayStlMultiset {
protected:
	struct SetDatum {
		float key;
		void *circularArrayDatum;
	};
	struct SetPredicate {
		Predicate predicate;
		bool operator () (const SetDatum& x,const SetDatum& y) const {
			return predicate(x.key,y.key);
		}
	};
	typedef std::multiset<SetDatum,SetPredicate> Set;
	struct Datum {
		typename Set::iterator iter;
	};
	CircularArray<Datum> circularArray;
	Set set;
	Predicate predicate;

	int64_t totalFullness;
	int totalEnters;
public:
	OpaqueQueueWithFindMin_CircularArrayStlMultiset (int maxSize) 
		: circularArray(maxSize)
	{
		totalFullness=0;
		totalEnters=0;
	}
	void resize (int maxSize) {
		circularArray.resize(maxSize);
	}
	void Enter (const T& key) {
		// remove lower-priority keys

		while (!set.empty()) {
			typename Set::iterator i=set.end();
			i--;
			if (!predicate(i->key,key)) {
				// i->key will be hidden for the rest of its life (in the case of equality, we only need one)
				Datum *datum=(Datum *)i->circularArrayDatum;
				datum->iter=set.end();
				set.erase(i);
			}
			else {
				break;
			}
		}

		Datum& slot=circularArray.EnterNull();
		SetDatum setDatum;
		setDatum.key=key;
		setDatum.circularArrayDatum=&slot;
		slot.iter=set.insert(setDatum);

		totalFullness += set.size();
		totalEnters++;
	}
	void Leave (void) {
		Datum& slot=circularArray.LeaveNull();
		if (slot.iter!=set.end()) {
			set.erase(slot.iter);
		}
	}
	const T& FindMin (void) const {
		return set.begin()->key;
	}
	void ClearAll (void) {
		if (totalEnters>0) {
			//printf("Average fullness = %lg\n",(double)(totalFullness)/(double)(totalEnters));
		}
		circularArray.ClearAll();
		set.clear();
		totalFullness=0;
		totalEnters=0;
	}
};

// if this works, I'm so smart
template <class T,class Predicate=std::less<T> >
class OpaqueQueueWithFindMin_CleverBinaryCircularArray {
protected:

	struct QueueDatum {
		int orderedPosition; // -1 means invalid
	};
	CircularArray<QueueDatum> circularArray;
	static int GetInvalidOrderedPosition (void) {
		return -1;
	}

	struct OrderedDatum {
		T key;
		QueueDatum *queueDatum;
	};
	vector<OrderedDatum> orderedArray;
	int circularFirst,circularLast;
	int size;
	Predicate predicate;
public:
	OpaqueQueueWithFindMin_CleverBinaryCircularArray (int maxSize) 
		: circularArray(maxSize)
	{
		size=maxSize;
		orderedArray.resize(maxSize);
		ClearAll();
	}
	void ClearAll (void) {
		circularArray.ClearAll();
		circularFirst=circularLast=0;
	}
	void resize (int maxSize) {
		size=maxSize;
		circularArray.resize(maxSize);
		orderedArray.resize(maxSize);
	}
	void Enter (const T& key) {

		// do a binary search in the orderedArray to find the position of our entry
		int insertPosition;
		if (circularFirst<=circularLast) {
			// simple binary search
			int first=circularFirst,last=circularLast;
			while (last-first>0) {
				int middle=(first+last)/2;
				if (predicate(orderedArray[middle].key,key)) {
					first=middle+1;
				}
				else {
					last=middle;
				}
			}
			insertPosition=last;
		}
		else {
			// weird wrapped binary search
			int first=circularFirst,last=circularLast+size;
			while (last-first>0) {
				int middle=(first+last)/2;
				int modMiddle=middle;
				if (modMiddle>=size) {
					modMiddle -= size;
				}
				if (predicate(orderedArray[modMiddle].key,key)) {
					first=middle+1;
				}
				else {
					last=middle;
				}
			}
			if (last>=size) {
				last -= size;
			}
			insertPosition=last;
		}

		// now, everything from insertPosition to circularLast is hidden, so we can remove it from orderedArray
		int i=insertPosition;
		while (i!=circularLast) {
			orderedArray[i].queueDatum->orderedPosition=GetInvalidOrderedPosition();
			i++;
			if (i==size) {
				i=0;
			}
		}
		circularLast=insertPosition;

		// and insert the new item
		QueueDatum& queueSlot=circularArray.EnterNull();
		orderedArray[circularLast].key=key;
		orderedArray[circularLast].queueDatum=&queueSlot;
		queueSlot.orderedPosition=insertPosition;
		
		circularLast++;
		if (circularLast==size) {
			circularLast=0;
		}
	}
	void Leave (void) {
		QueueDatum& slot=circularArray.LeaveNull();
		if (slot.orderedPosition!=GetInvalidOrderedPosition()) {

			// This is so cool!  I didn't realize it, but it's guaranteed that the item that's leaving the queue is the highest-priority item, which means that it's guaranteed to be in position circularFirst in the ordered array.  Why must the now-leaving item have highest priority?  Well, if there was some other element in the queue that had higher priority, it would have hidden the now-leaving item, so the now-leaving-item would be marked invalid
			assert(slot.orderedPosition==circularFirst);

			circularFirst++;
			if (circularFirst==size) {
				circularFirst=0;
			}

			/* Old code before I realized this wonderfulness
			int newPosition=slot.orderedPosition;
			int oldPosition=newPosition+1;
			if (oldPosition==size) {
				oldPosition=0;
			}
			while (oldPosition!=circularLast) {
				QueueDatum& queueSlot=*(orderedArray[oldPosition].queueDatum);
				queueSlot.orderedPosition=newPosition;
				orderedArray[newPosition]=orderedArray[oldPosition];
				newPosition=oldPosition;
				oldPosition++;
				if (oldPosition==size) {
					oldPosition=0;
				}
			}
			circularLast--;
			if (circularLast<0) {
				circularLast=size-1;
			}
			*/
		}
	}
	const T& FindMin (void) const {
		return orderedArray[circularFirst].key;
	}
};
