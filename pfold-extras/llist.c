/*********************************************************************

  Contains general functions for handling linked lists.

  990112 Bjarne Knudsen (bk@daimi.au.dk)

  Copyright (C) 1999 Bjarne Knudsen

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

*********************************************************************/

#include <stdlib.h>
#include "llist.h"

/* Allocate and initialize LList */
LList *MakeLList()
{
  LList *llist = (LList *)malloc(sizeof(LList));

  /* Initialize llist */
  llist->first = (LListNode *)malloc(sizeof(LListNode));
  llist->last = (LListNode *)malloc(sizeof(LListNode));
  llist->count = 0;

  /* initialize stop blocks */
  (llist->first)->elm = NULL;
  (llist->first)->prev = NULL;
  (llist->first)->next = llist->last;
  (llist->last)->elm = NULL;
  (llist->last)->prev = llist->first;
  (llist->last)->next = NULL;

  return llist;
}

void DestroyLList(LList *llist)
{
  LListNode *lnode;

  /* free LListNodes */
  while ((lnode = llist->first->next) != NULL) {
    (llist->first)->next = lnode->next;
    free(lnode);
  }
  free(llist->first);

  /* free LList */
  free(llist);
}

/* Insert node after 'lnode' */
void Insert(LList *llist, LListNode *lnode, void *elm)
{
  LListNode *newnode = (LListNode *)malloc(sizeof(LListNode));

  /* initialize newnode */
  newnode->next = lnode->next;
  newnode->prev = lnode;
  newnode->elm = elm;

  /* insert it */
  (newnode->prev)->next = newnode;
  (newnode->next)->prev = newnode;
  (llist->count)++;
}

/* Remove 'lnode' */
void *Remove(LList *llist, LListNode *lnode)
{
  void *elm = lnode->elm;

  (lnode->prev)->next = lnode->next;
  (lnode->next)->prev = lnode->prev;
  (llist->count)--;
  free(lnode);

  return elm;
}

void Push(LList *llist, void *elm)
{
  Insert(llist, llist->first, elm);
}

void *Pop(LList *llist)
{
  void *elm = NULL;

  if ((llist->first)->next != llist->last)
    /* Not empty list */
    elm = Remove(llist, (llist->first)->next);

  return elm;
}

void Enqueue(LList *llist, void *elm)
{
  Insert(llist, (llist->last)->prev, elm);
}

void *Dequeue(LList *llist)
{
  void *elm = NULL;

  if ((llist->last)->prev != llist->first)
    /* Not empty list */
    elm = Remove(llist, (llist->last)->prev);

  return elm;
}

void LListMap(LList *llist, void (*f)(void *, void *), void *arg)
{
  LListNode *lnode;

  for(lnode = (llist->first)->next; lnode->next != NULL; lnode = lnode->next)
    (*f)(lnode->elm, arg);
}

LListCounter *MakeCounter(LList *llist, int pos)
{
  LListCounter *lcounter = (LListCounter *)malloc(sizeof(LListCounter));

  InitCounter(lcounter, llist, pos);

  return lcounter;
}

void InitCounter(LListCounter *lcounter, LList *llist, int pos)
{
  if (pos == LAST)
    pos = llist->count;

  if (pos < llist->count - 1 - pos) {
    /* closest to first node */
    lcounter->pos = 0;
    lcounter->current = (llist->first)->next;
  }
  else {
    /* closest to last node */
    lcounter->pos = llist->count - 1;
    lcounter->current = (llist->last)->prev;
  }
  SetCounter(lcounter, pos);
}

void *SetCounter(LListCounter *lcounter, int pos)
{
  if(lcounter->pos < pos)
    /* move lcounter forward */
    while(lcounter->pos < pos && (lcounter->current)->next != NULL) {
      (lcounter->pos)++;
      lcounter->current = (lcounter->current)->next;
    }
  else if (lcounter->pos > pos)
    /* move lcounter backward */
    while(lcounter->pos > pos && (lcounter->current)->prev != NULL) {
      (lcounter->pos)--;
      lcounter->current = (lcounter->current)->prev;
    }

  return (lcounter->current)->elm;
}

void *Next(LListCounter *lcounter)
{
  return SetCounter(lcounter, lcounter->pos + 1);
}

void *Prev(LListCounter *lcounter)
{
  return SetCounter(lcounter, lcounter->pos - 1);
}
