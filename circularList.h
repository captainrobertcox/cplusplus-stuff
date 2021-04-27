// Robert Cox, B00226347, cox.robert4@titans.easternflorida.edu,
// COP2335-70B, Module 8, Creating a circular linked list, 4-17-21

#ifndef LinkedList_H
#define LinkedList_H

#include <iostream>
#include <string>

using namespace std;

template <class Type>
struct node{
  Type info;
  node<Type>* link;
};

template <class Type>
class linkedListIterator{

  public:

    linkedListIterator(); // default constructor
      // Post-condition: current = nullptr;

    linkedListIterator(node<Type>* ptr);
      // Constructor with parameter
      // Pre-condition: ptr must be a valid node<Type> pointer
      // Post-condition: current = ptr

    Type operator*(); // overload dereferencing operator
      // Post-condition: Returns info of the node

    linkedListIterator<Type> operator++(); // overload pre-increment
      // overloaded pre-increment operator
      // Post-condition: the iterator is advanced to the next node

    linkedListIterator<Type> operator++(int); // overload post-increment
      // overloaded post-increment operator
      // Post-condition: the iterator is advanced to the next node after this operation

    bool operator==(const linkedListIterator<Type>& right) const;
      // overloaded equality operator
      // Pre-condition: right iterator must be a valid linkedListIterator<Type>
      // Post-condition: returns true if this iterator is equal to the right iterator, false otherwise.

    bool operator!=(const linkedListIterator<Type>& right) const;
      // overloaded non-equality operator
      // Pre-condition: right iterator must be a valid linkedListIterator<Type>
      // Post-condition: returns true if this iterator is not equal to the right iterator, false if they are equal.

  private:

    node<Type>* current; // points to the current node in the linked list
};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
linkedListIterator<Type>::linkedListIterator(){
  current = nullptr;
}

//--------------------------------------------------------------------------------

template <class Type>
linkedListIterator<Type>::linkedListIterator(node<Type>* ptr){
  current = ptr;
}

//--------------------------------------------------------------------------------

template <class Type>
Type linkedListIterator<Type>::operator*(){
  return current->info;
}

//--------------------------------------------------------------------------------

template <class Type>
linkedListIterator<Type> linkedListIterator<Type>::operator++(){
  current = current->link;
  return *this;
}

//--------------------------------------------------------------------------------

template <class Type>
linkedListIterator<Type> linkedListIterator<Type>::operator++(int u){
  linkedListIterator<Type> temp = *this;
  current = current->link;
  return temp;
}

//--------------------------------------------------------------------------------

template <class Type>
bool linkedListIterator<Type>::operator==(const linkedListIterator& right) const{
  return (current == right.current);
}

//--------------------------------------------------------------------------------

template <class Type>
bool linkedListIterator<Type>::operator!=(const linkedListIterator& right) const{
  return (current != right.current);
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
class linkedList{

  public:

    const linkedList<Type>& operator=(const linkedList<Type>&);
      // overloaded assignment operator
      // Pre-condition: other list must be a valid linkedList<Type>
      // Post-condition: this list is copied from the other list

    bool operator==(const linkedList<Type>&) const;
      // overloaded equality operator
      // Pre-condition: other list must be a valid linkedList<Type>
      // Post-condition: returns true if this list has the same info as the other list, false otherwise

    void initializeList();
      // initialize the list to an empty state
      // Post-conditions: last = nullptr, count = 0

    bool isEmpty() const;
      // determines if the list is empty
      // Post-condition: returns true if the list is empty, false otherwise.

    void print() const;
      // this function displays the info in each node of the list
      // Post-condition: prints the list

    int length() const;
      // this function returns the number of nodes in the list
      // Post-condition: returns the value of count, which is the number of nodes in the list.

    void destroy();
      // this function deletes all the nodes in the list
      // Post-conditions: last = nullptr, count = 0, all the nodes are deleted.

    Type front() const;
      // this function returns the first element in the list
      // Pre-condition: the list must exist and not be empty.
      // Post-condition: if the list is empty, a string exception is thrown. Otherwise, the first node is returned.

    Type back() const;
      // this function returns the last element in the list
      // Pre-condition: the list must exist and not be empty
      // Post-condition: if the list is empty, a string exception is thrown. Otherwise, the last node is returned.

    virtual bool search(const Type& searchItem) const = 0; // pure virtual with runtime binding
      // function to determine if searchItem is in the list
      // Pre-condition: searchItem must be a valid Type for the list
      // Post-condition: returns true if searchItem is in the list, false otherwise.

    virtual void insertFirst(const Type& newItem) = 0;
      // function to insert newItem at beginning of the list
      // Pre-condition: newItem must be a valid Type for the list
      // Post-condition: newItem is inserted at the beginning of the list, count is incremented.

    virtual void insertLast(const Type& newItem) = 0;
      // function to insert newItem at end of the list
      // Pre-condition: newItem must be a valid Type for the list
      // Post-condition: newItem is inserted at the end of the list, count is incremented.

    virtual void deleteNode(const Type& deleteItem) = 0;
      // function to delete a node from the list
      // Pre-condition: deleteItem must be a valid Type for the list.
      // Post-condition: if deleteItem is in the list, it is removed from the list, and count is decremented.

    linkedListIterator<Type> begin() const;
      // this function returns an iterator for the beginning of the list
      // Post-condition: returns an iterator, whose current is set at the beginning of the list

    //linkedListIterator<Type> end() const;
      // This does not have meaning for a circular list, so it is not used.
      // this function returns an iterator for the end of the list. In a circular list, the end is the beginning.
      // Post-condition: returns an iterator, whose current is set to the end of the list.

    linkedList(); // constructor
      // default constructor
      // Post-conditions: initializes list to an empty state, sets count = 0, last = nullptr.

    linkedList(const linkedList<Type>& otherList); // overload copy constructor
      // constructor with otherList as a parameter
      // Pre-condition: otherList must be a valid linkedList<Type>
      // Post-conditions: copies otherList info into this list

    ~linkedList(); // destructor
      // destructor
      // Post-conditions: destroys the list by deleting all the nodes. count will be 0.

  protected:

    int count; // holds size of list
    //node<Type>* first; // for a circular list, last->link is first
    node<Type>* last;

  private:

    void copyList(const linkedList<Type>& otherList);
};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
bool linkedList<Type>::isEmpty() const {
  return (last == nullptr);
}

//--------------------------------------------------------------------------------

template <class Type>
bool linkedList<Type>::operator==(const linkedList<Type>& otherList) const{
  if (otherList.count != count) return false;
  node<Type>* current = last;
  node<Type>* otherListCurrent = otherList.last;
  for (int i = 0; i < count; i++){
    if (current->info != otherListCurrent->info) return false;
    current = current->link;
    otherListCurrent = otherListCurrent->link;
  }
  return true;
}

//--------------------------------------------------------------------------------

template <class Type>
linkedList<Type>::linkedList(){
  last = nullptr;
  count = 0;
}

//--------------------------------------------------------------------------------

template <class Type>
void linkedList<Type>::destroy(){
  node<Type>* temp;
  for (int i = 0; i < count; i++){
    temp = last;
    last = last->link;
    delete temp;
    temp = nullptr;
  }
  last = nullptr;
  count = 0;
}

//--------------------------------------------------------------------------------

template <class Type>
void linkedList<Type>::initializeList(){
  destroy();
}

//--------------------------------------------------------------------------------

template <class Type>
void linkedList<Type>::print() const{
  if (count == 0) return;
  node<Type>* current;
  current = last->link; // first;
  for (int i = 0; i < count; i++){
    cout << current->info << " ";
    current = current->link;
  }
}

//--------------------------------------------------------------------------------

template <class Type>
int linkedList<Type>::length() const {
  return count;
}

//--------------------------------------------------------------------------------

template <class Type>
Type linkedList<Type>::front() const{
  try{
    if (last == nullptr) throw string("List is empty");
    if (last->link != nullptr) return last->link->info; // error if there is no first node yet
  }
  catch (string str){
    throw str;
  }
}

//--------------------------------------------------------------------------------

template <class Type>
Type linkedList<Type>::back() const{
  try{
    if (last == nullptr) throw string("List is empty");
    if (last != nullptr) return last->info; // before calling this function, be certain not empty
  }
  catch (string str){
    throw str;
  }
}

//--------------------------------------------------------------------------------

template <class Type>
linkedListIterator<Type> linkedList<Type>::begin() const {
  linkedListIterator<Type> temp(last->link); // first
  return temp;
}

//--------------------------------------------------------------------------------

/*
template <class Type>
linkedListIterator<Type> linkedList<Type>::end() const { // pointless function - same as begin
  linkedListIterator<Type> temp(last->link); // is there and end in a circular list
  return temp;
}
*/

//--------------------------------------------------------------------------------

template <class Type>
void linkedList<Type>::copyList(const linkedList<Type>& otherList){
  node<Type>* newNode;
  node<Type>* current;
  if (last != nullptr) destroy();
  if (otherList.last == nullptr){
    last = nullptr;
    count = 0;
  }
  else{
    current = (otherList.last)->link; // first;
    count = otherList.count;
    last = new node<Type>;
    last->info = current->info;
    last->link = last;
    current = current->link;
    for (int i = 0; i < count-1; i++){
      newNode = new node<Type>;
      newNode->info = current->info;
      newNode->link = last->link;
      last->link = newNode;
      last = newNode;
      current = current->link;
    } // end for
  } // end else
}

//--------------------------------------------------------------------------------

template <class Type>
linkedList<Type>::~linkedList(){
  destroy();
}

//--------------------------------------------------------------------------------

template <class Type>
linkedList<Type>::linkedList(const linkedList<Type>& otherList){
  last = nullptr;
  copyList(otherList);
}

//--------------------------------------------------------------------------------

template <class Type>  // overload assignment operator
const linkedList<Type>& linkedList<Type>::operator=(const linkedList<Type>& otherList){
  if (this != &otherList){ // no self assignment
    copyList(otherList);
  }
  return *this;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
class unorderedList : public linkedList<Type> {

  public:

    bool search(const Type& searchItem) const;
      // this function determines if searchItem is in the list
      // Pre-conditions: searchItem must be a valid Type for the list
      // Post-conditions: returns true if searchItem is in the list, false otherwise.

    void insertFirst(const Type& newItem);
      // function to insert newItem at beginning of the list
      // Pre-conditions: newItem must be a valid Type for the list
      // Post-conditions: inserts newItem at beginning of the list, increments count, last->link is newItem.

    void insertLast(const Type& newItem);
      // function to insert newItem at end of list
      // Pre-conditions: newItem must be a valid Type for the list
      // Post-conditions: newItem is inserted at the end of the list, increments count, last is newItem.

    void deleteNode(const Type& deleteItem);
      // function to delete a node from the list
      // Pre-condition: deleteItem must be a valid Type for the list
      // Post-condition: if deleteItem is in the list, it is removed, and count is decremented.

};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
bool unorderedList<Type>::search(const Type& searchItem) const{
  node<Type>* current;
  current = this->last;
  for (int i = 0; i < this->count; i++){
    if (current->info == searchItem) return true;
    else current = current->link;
  }
  return false;
}

//--------------------------------------------------------------------------------

template <class Type>
void unorderedList<Type>::insertFirst(const Type& newItem){
  node<Type>* newNode;
  newNode = new node<Type>;
  newNode->info = newItem;
  if (this->last == nullptr){
    this->last = newNode;
    this->last->link = newNode;
    this->count++;
  }
  else{
    newNode->link = this->last->link;
    this->last->link = newNode;
    this->count++;
  }
}

//--------------------------------------------------------------------------------

template <class Type>
void unorderedList<Type>::insertLast(const Type& newItem){
  node<Type>* newNode;
  newNode = new node<Type>;
  newNode->info = newItem;
  newNode->link = this->last->link;
  if (this->last == nullptr){
    this->last = newNode;
    this->last->link = newNode;
    this->count++;
  }
  else{
    newNode->link = this->last->link;
    this->last->link = newNode;
    this->last = newNode;
    this->count++;
  }
}

//--------------------------------------------------------------------------------

template <class Type>
void unorderedList<Type>::deleteNode(const Type& deleteItem){
  node<Type>* current;
  node<Type>* trailCurrent;
  bool found;
  if (this->last == nullptr) cout << "Cannot delete from empty list." << endl;
  else{
    if (this->last->link->info == deleteItem){
      current = this->last->link;
      this->last->link = this->last->link->link;
      this->count--;
      if (this->count == 0) this->last = nullptr;
      delete current;
    }
    else{
      found = false;
      trailCurrent = this->last;
      current = this->last->link;
      for (int i = 0; i < this->count && !found; i++){
        if (current->info != deleteItem){
          trailCurrent = current;
          current = current->link;
        }
        else found = true;
      } // for
      if (found){
        trailCurrent->link = current->link;
        this->count--;
        if (this->last == current) this->last = trailCurrent;
        delete current;
      }
      else cout << "The item is not in the list." << endl;
    } // else
  } // else
}

//-----------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
class orderedList : public linkedList<Type>{

  public:


    bool search(const Type& searchItem) const;
      // function to determine if searchItem is in the list
      // Pre-condition: searchItem must be a valid Type for the list
      // Post-condition: returns true if searchItem is in the list, false otherwise.

    void insert(const Type& newItem);
      // function to insert newItem in the correct place in the list
      // Pre-condition: newItem must be a valid Type for the list
      // Post-condition: newItem is inserted into the correct position in the list, count is incremented.

    void insertFirst(const Type& newItem);
      // for ordered lists, insert function is called to insert newItem in correct place in list
      // Pre-condition: newItem must be a valid Type for the list
      // Post-condition: insert(newItem) is called, which inserts item in correct place and increments count.

    void insertLast(const Type& newItem);
      // for ordered lists, insert function is called to insert newItem in correct place in list
      // Pre-condition: newItem must be a valid Type for the list
      // Post-condition: insert(newItem) is called, which inserts item in correct place and increments count.

    void deleteNode(const Type& deleteItem);
      // function to remove deleteItem from the list
      // Pre-condition: deleteItem must be a valid Type for the list
      // Post-conditions: if deleteItem is in the list, it is removed and count is decremented.

};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <class Type>
bool orderedList<Type>::search(const Type& searchItem) const{
  bool found = false;
  node<Type>* current;
  current = this->last->link;
  for (int i = 0; i < this->count && !found; i++){
    if (current->info >= searchItem) found = true;
    else current = current->link;
  }
  if (found) found = (current->info == searchItem);
  return found;
}

//--------------------------------------------------------------------------------

template <class Type>
void orderedList<Type>::insert(const Type& newItem){
  node<Type>* current;
  node<Type>* trailCurrent = nullptr;
  node<Type>* newNode;
  bool found;
  newNode = new node<Type>;
  newNode->info = newItem;
  if (this->last == nullptr){
    this->last = newNode;
    this->last->link = newNode;
    this->count++;
  }
  else{
    current = this->last->link;
    trailCurrent = this->last;
    found = false;
    for (int i = 0; i < this->count && !found; i++){
      if (current->info >= newItem) found = true;
      else{
        trailCurrent = current;
        current = current->link;
      }
    } // for
    if (!found){
      // new item goes at end of list
      newNode->link = this->last->link;
      trailCurrent->link = newNode;
      this->last = newNode;
      this->count++;
    }
    else{
    if (current == this->last){
      // item is inserted just before the last
      newNode->link = this->last;
      trailCurrent->link = newNode;
      this->count++;
    }
    else{
      if (current == this->last->link) this->last->link = newNode; // newNode is new first
      else{
        trailCurrent->link = newNode;
      }
      newNode->link = current;
      this->count++;
    }
    }
  } // else
}

//--------------------------------------------------------------------------------

template <class Type>
void orderedList<Type>::insertFirst(const Type& newItem){
  insert(newItem);
}

//--------------------------------------------------------------------------------

template <class Type>
void orderedList<Type>::insertLast(const Type& newItem){
  insert(newItem);
}

//--------------------------------------------------------------------------------

template <class Type>
void orderedList<Type>::deleteNode(const Type& deleteItem){
  node<Type>* current;
  node<Type>* trailCurrent;
  bool found;
  if (this->last == nullptr) cout << "Cannot delete from empty list." << endl;
  else{
    current = this->last->link;
    trailCurrent = this->last;
    found = false;
    for (int i = 0; i < this->count && !found; i++){
      if (current->info >= deleteItem) found = true;
      else{
        trailCurrent = current;
        current = current->link;
      }
    }
    if (!found) cout << "The item is not in the list," << endl;
    else{
      if (current->info == deleteItem){
        if (this->last->link == current){
          // delete first item in list
          if (this->count == 1) this->last = nullptr;
          else{
            this->last->link = this->last->link->link;
          }
          delete current;
        }
        else{
          trailCurrent->link = current->link;
          if (current == this->last) this->last = trailCurrent;
          delete current;
        }
        this->count--;
      } // if
      else cout << "The item is not in the list." << endl;
    } // else
  } // else
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

#endif
