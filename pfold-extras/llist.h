/*********************************************************************
 *                                                                   *
 * This is an interface for a linked list as described in llist.txt. *
 * This draws heavily from work by C. Storm.                         *
 *                                                                   *
 * 990112 Bjarne Knudsen                                             *
 *                                                                   *
 *********************************************************************/

#ifndef _LLIST_H
#define _LLIST_H

#define FIRST  -1
#define LAST   -2

typedef struct tagLListNode {
  void *elm;
  struct tagLListNode *prev;
  struct tagLListNode *next;
} LListNode;

typedef struct tagLList {
  int count;
  struct tagLListNode *first;
  struct tagLListNode *last;
} LList;

typedef struct tagLListCounter {
  int pos;
  struct tagLListNode *current;
} LListCounter;

/* prototypes */

LList *MakeLList(void);
void DestroyLList(LList *llist);
void Insert(LList *llist, LListNode *lnode, void *elm);
void *Remove(LList *llist, LListNode *lnode);
void Push(LList *llist, void *elm);
void *Pop(LList *llist);
void Enqueue(LList *llist, void *elm);
void *Dequeue(LList *llist);
void LListMap(LList *llist, void (*f)(void *, void *), void *arg);
LListCounter *MakeCounter(LList *llist, int);
void InitCounter(LListCounter *lcounter, LList *llist, int pos);
void *SetCounter(LListCounter *lcounter, int pos);
void *Next(LListCounter *lcounter);
void *Prev(LListCounter *lcounter);

#endif
