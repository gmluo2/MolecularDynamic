template <class T> class LIST {
public:
	T *p;
	LIST<T> *next, *pre;
	LIST() {p = NULL; next = NULL; pre = NULL;};
	~LIST() {p = NULL; next = NULL; pre = NULL;};
	void Attach(LIST<T> *l) {next = l; l->pre = this;};
	LIST<T>* get_list(T* t) {
		LIST<T> *list = this;
		while (list != NULL && list->p != t) list = list->next;
		return list;
	};
	LIST<T>* get_head() {
		LIST<T> *p = this;
		while (p->pre != NULL) p = p->pre;
		return p;
	};
	LIST<T>* get_tail() {
		LIST<T> *p = this;
		while (p->next != NULL) p = p->next;
		return p;
	};
	LIST<T>* get_list(int n) {
		LIST<T> *p = NULL;
		int i = 0;
		if (n >= 0) {
			p = this->get_head();
			while (i < n && p != NULL) {p = p->next; i++;}
			return p;
		}
		else {
			n = -n;
			p = this->get_tail();
			while (i < n) {p = p->pre; i++;}
			return p;
		}
	};
	int nList(bool empty = false) {
		LIST<T> *plist = this;
		int n = 0;
		while (plist != NULL) {
			if (empty) n++;
			else if (plist->p != NULL) n++;
			plist = plist->next;
		}
		return n;
	};
	int indx(LIST<T> *p) {
		LIST<T> *plist = this;
		int n = 0;
		while (plist != p) {
			n++; plist = plist->next;
			if (plist == NULL) return -1;
		}
		return n;
	};
	int indx(T *p) {
		LIST<T> *plist = this;
		int n = 0;
		while (plist->p != p) {
			n++; plist = plist->next;
			if (plist == NULL) return -1;
		}
		return n;
	};
	LIST<T>* add(T *t, bool duplicatable) {
		LIST<T> *plist = this->get_head();
		LIST<T> *p = NULL;
		if (!duplicatable && (p = plist->get_list(t)) != NULL) return p;
		while (plist->next != NULL && plist->p != NULL) plist = plist->next;
		if (plist->p == NULL) {plist->p = t; return plist;}
		else { // plist->p != NULL, so plist must be the tail
			p = new LIST<T>; p->p = t;
			plist->Attach(p);
			return p;
		}
	};
};

template <class T> void release_list(LIST<T> **list, bool release_context) {
	if (*list == NULL) return;
	LIST<T> *next = NULL, *plist = *list;
	if (plist->pre != NULL) plist->pre->next = NULL;
	while (plist != NULL) {
		next = plist->next;
		if (release_context && plist->p != NULL) {delete plist->p; plist->p = NULL;}
		else plist->p = NULL;
		delete plist;
		plist = next;
	}
	*list = NULL;
};

template <class T> void detach(LIST<T> **list, LIST<T> *pdetach, bool release_context) {
	LIST<T> *pre = pdetach->pre, *next = pdetach->next;
	if (pdetach == *list) *list = pdetach->next;
	if (pre != NULL) pre->next = next;

	if (next != NULL) next->pre = pre;
	if (release_context) delete pdetach->p;
	pdetach->p = NULL;
	delete pdetach;
}

template <class T> void detach(LIST<T> **list, T *pdetach_T, bool release_context) {
	LIST<T> *pdetach = (*list)->get_list(pdetach_T);
	LIST<T> *pre = pdetach->pre, *next = pdetach->next;
	if (pdetach == *list) *list = pdetach->next;
	if (pre != NULL) pre->next = next;

	if (next != NULL) next->pre = pre;
	if (release_context) delete pdetach->p;
	pdetach->p = NULL;
	delete pdetach;
}

template <class T> class TREE {
public:
	LIST< LIST<T> > *tl; // the level from the root
	int nlevels;

	TREE() {tl = NULL; nlevels = 0;};
	~TREE() {tl = NULL; nlevels = 0;};
	void release_tree(bool release_context) {
		LIST< LIST<T> > *ptl = tl;
		while (ptl != NULL) {
			release_list<T>(&(ptl->p), release_context);
			ptl = ptl->next;
		}
		release_list< LIST<T> >(&tl, false); // the context inside the tree levels are clean
		nlevels = 0;
	};
	LIST<T>* get_level(int n) {
		if (n < 0) n += nlevels; // start from the last level
		LIST< LIST<T> > *ptl = tl;
		int i = 0;
		while (i < n && ptl != NULL) {ptl = ptl->next; i++;};
		if (ptl == NULL) return ptl;
		else return ptl->p;
	};
	LIST< LIST<T> >* get_tree_level(int n) {
		if (n < 0) n += nlevels; // start from the last level
		LIST< LIST<T> > *ptl = tl;
		int i = 0;
		while (i < n && ptl != NULL) {ptl = ptl->next; i++;};
		if (ptl == NULL) return ptl;
		else return ptl;
	};
	LIST< LIST<T> >* new_level() {
		LIST< LIST<T> > *p = new LIST< LIST<T> >;
		if (tl == NULL) {tl = p; return p;}
		LIST< LIST<T> > *tail = tl->get_tail();
		tail->Attach(p);
		nlevels++;
		return p;
	};
	LIST<T>* add(LIST< LIST<T> >* ptl, T* t, bool duplicatable) {
		LIST<T> *p = NULL;
		if (ptl->p == NULL) {ptl->p = new LIST<T>; ptl->p->p = t; return p;}
		p = ptl->p;
		return p->add(t, duplicatable);
	};
};

template <class T> class CHAIN {
public:
	T* p;
	CHAIN *next;
	CHAIN() {p = NULL; next = NULL;};
	~CHAIN() {p = NULL; next = NULL;};
	void attach2tail(T* p, bool non_repeat = true) {
		CHAIN<T> *tch = NULL, *ch = NULL;
		ch = this; 
		if (non_repeat && ch->p == p) return; // p is in the chain
		while (ch->next != NULL) {
			ch = ch->next;
			if (non_repeat && ch->p == p) return; // p is in the chain
		}
		tch = new CHAIN<T>; tch->p = p;
		ch->next = tch;
	};
	void insert(T* p) {
		CHAIN<T> *tch = new CHAIN<T>, *ch = this->next;
		this->next = tch; tch->next = ch;
	};
	CHAIN<T>* get_tail() {
		CHAIN<T> *ch = this;
		while (ch->next != NULL) ch = ch->next;
		return ch;
	};
};

template <class T> void release_chain(CHAIN<T> **ch, bool release_context = false) {
	CHAIN<T> *pch = *ch, *next = NULL;
	while (pch != NULL) {
		next = pch->next; pch->next = NULL;
		if (release_context) {delete pch->p; pch->p = NULL;}
		delete pch;
		pch = next;
	}
	*ch = NULL;
};

template <class T> bool in_chain(CHAIN<T> *ch, T* p) {
	CHAIN<T> *pch = ch;
	while (pch != NULL) {
		if (pch->p == p) return true;
		pch = pch->next;
	}
	return false;
};

template <class T> bool exist_in_chain(CHAIN<T> *ch, T& p) { // context is the same
	CHAIN<T> *pch = ch;
	while (pch != NULL) {
		if (*(pch->p) == p) return true;
		pch = pch->next;
	}
	return false;
};

template <class T> bool detach_from_chain(CHAIN<T> **ch, T* t, bool release_context) {
	CHAIN<T> *pch = *ch, *parent = *ch;
	while (pch != NULL && pch->p != t) {
		parent = pch;
		pch = pch->next;
	}
	if (pch != NULL) {
		if (pch == *ch) *ch = pch->next;
		else parent->next = pch->next;
		if (release_context) {delete pch->p; pch->p = NULL;}
		delete pch;
		return true;
	}
	else return false;
};

template <class T> int number(CHAIN<T> *ch) {
	int n = 0;
	CHAIN<T> *p = ch;
	while (p != NULL) {n++; p = p->next;}
	return n;
};

template <class T> bool in_tree(T *atom, TREE<T> &atree) {
	LIST< LIST<T> > *plevel = atree.tl;
	LIST<T> *plist = NULL;

	while (plevel != NULL) {
		plist = plevel->p;
		while (plist != NULL && plist->p != NULL) {
			if (plist->p == atom) return true;
			plist = plist->next;
		}
		plevel = plevel->next;
	}
	return false;
}

template <class T> bool search(T *atom, TREE<T> &atree, LIST< LIST<T> > **level, LIST<T> **list) {
	LIST< LIST<T> > *plevel = atree.tl;
	LIST<T> *plist = NULL;

	*level = NULL; *list = NULL;

	while (plevel != NULL) {
		plist = plevel->p;
		while (plist != NULL && plist->p != NULL) {
			if (plist->p == atom) {*level = plevel; *list = plist; return true;}
			plist = plist->next;
		}
		plevel = plevel->next;
	}
	return false;
}

/*************************************************************
   class WRAP has an element T* pbase
   try to find out the WRAP list has element p inside
*************************************************************/
template <class WRAP, class T> LIST<WRAP>* get_wrap_list(LIST<WRAP> *plist, T* p) {
	LIST<WRAP> *pw = plist->get_head();
	while (pw != NULL && pw->p->pbase != p) pw = pw->next;
	return pw;
};

#include "dstruct1.h"

template <class T> void array_chain(ARRAY<T*> &a, CHAIN<T> **ch) {
	if (ch != NULL) {
		if (*ch != NULL) release_chain<T>(ch, true);
	}
	if (a.n == 0 || a.m == NULL) return;
	CHAIN<T> *tail = NULL, *pch = NULL;
	for (int i = 0; i < a.n; i++) {
		pch = new CHAIN<T>; pch->p = a.m[i];
		if (*ch == NULL) *ch = pch;
		else tail->next = pch;
		tail = pch;
	}
};

template <class T> void chain_array(CHAIN<T> **ch, ARRAY<T*> &a) {
	a.SetArray(0);
	if (ch == NULL || *ch == NULL) return;
	int n = number<T>(*ch), i = 0;
	if (n == 0) return;
	a.SetArray(n);
	CHAIN<T> *pch = *ch;
	while (pch != NULL) {
		a.m[i] = pch->p;
		pch = pch->next; i++;
	}
};

template <class T> void chain_array_context_cp(CHAIN<T> **ch, ARRAY<T> &a) {
	a.SetArray(0);
	if (ch == NULL || *ch == NULL) return;
	int n = number<T>(*ch), i = 0;
	if (n == 0) return;
	a.SetArray(n);
	CHAIN<T> *pch = *ch;
	while (pch != NULL) {
		memcpy(a.m + i, pch->p, sizeof(T)); //a.m[i] = pch->p;
		pch = pch->next; i++;
	}
};
