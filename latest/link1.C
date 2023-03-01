// FILE: link1.cxx
// IMPLEMENTS: The ten functions of the linked list toolkit (see link1.h)

#include <assert.h>    // Provides assert
#include <stdlib.h>    // Provides NULL and size_t
#include "link1.h"

void list_head_insert(Node*& head_ptr, const Node::Item& entry)
{
    Node *insert_ptr;

    insert_ptr = new Node;
    insert_ptr->data = entry;
    insert_ptr->link = head_ptr;
    head_ptr = insert_ptr;
}

void list_insert(Node* previous_ptr, const Node::Item& entry) 
{
    Node *insert_ptr;
    
    insert_ptr = new Node;
    insert_ptr->data = entry;
    insert_ptr->link = previous_ptr->link;
    previous_ptr->link = insert_ptr;
}

void list_head_remove(Node*& head_ptr)
{
    Node *remove_ptr;

    remove_ptr = head_ptr;
    head_ptr = head_ptr->link;
    delete remove_ptr;
}

void list_remove(Node* previous_ptr)
{
    Node *remove_ptr;

    remove_ptr = previous_ptr->link;
    previous_ptr->link = remove_ptr->link;
    delete remove_ptr;
}

void list_clear(Node*& head_ptr)
// Library facilities used: stdlib.h
{
    while (head_ptr != NULL)
        list_head_remove(head_ptr);
}

void list_copy(Node* source_ptr, Node*& head_ptr, Node*& tail_ptr) 
// Library facilities used: stdlib.h
{
    head_ptr = NULL;
    tail_ptr = NULL;

    // Handle the case of the empty list
    if (source_ptr == NULL)
        return;
    
    // Make the head node for the newly created list, and put data in it
    list_head_insert(head_ptr, source_ptr->data);
    tail_ptr = head_ptr;
    
    // Copy the rest of the nodes one at a time, adding at the tail of new list
    for (source_ptr = source_ptr->link; source_ptr != NULL; source_ptr = source_ptr->link)
    {
        list_insert(tail_ptr, source_ptr->data);
        tail_ptr = tail_ptr->link;
    }
}

int list_length(Node *head_ptr)
{
  int length = 0;
  for(;head_ptr!=NULL;head_ptr=head_ptr->link)
    length++;
  return length;
}
