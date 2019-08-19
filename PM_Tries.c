/*trie.c*/
#include "PM.h"


void TrieCreate(trieNode_t **root)
{
   *root = TrieCreateNode('\0', 0xffffffff);
}

trieNode_t *TrieCreateNode(char key, int data)
{
   trieNode_t *node = NULL;
   node = (trieNode_t *)malloc(sizeof(trieNode_t));

   if (NULL == node)
   {
      printf("Malloc failed\n");
      return node;
   }

   node->key = key;
   node->next = NULL;
   node->children = NULL;
   node->value = data;
   node->count = NULL;
   node->parent = NULL;
   node->prev = NULL;
   node->prev = NULL;
   return node;
}

trieValueNode_t *TrieCreateValueNode(int data)
{
   trieValueNode_t *node = NULL;
   node = (trieNode_t *)malloc(sizeof(trieValueNode_t));

   if (NULL == node)
   {
      printf("Malloc failed\n");
      return node;
   }
   node->value = data;
   return node;
}


void TrieAdd(trieNode_t **root, char *key, int data, posList_t **headpos)
{
   trieNode_t *pTrav = NULL;

   if (NULL == *root)
   {
      printf("NULL tree\n");
      return;
   }
#ifdef DEBUG
   printf("\nInserting key %s: \n", key);
#endif
   pTrav = (*root)->children;

   if (pTrav == NULL)
   {
      /*First Node*/
      for (pTrav = *root; *key; pTrav = pTrav->children)
      {
         pTrav->children = TrieCreateNode(*key, 0xffffffff);
         pTrav->children->parent = pTrav;
#ifdef DEBUG
         printf("Inserting: [%c]\n", pTrav->children->key);
#endif
         key++;
      }

      pTrav->children = TrieCreateNode('\0', 0);
      pTrav->children->parent = pTrav;
#ifdef DEBUG
      printf("Inserting: [%c]\n", pTrav->children->key);
#endif
      return;
   }

   trieNode_t *pPtr = TrieSearch(pTrav, key);
   if (pPtr)         //pPtr points at last node
   {
      //printf("Duplicate!\n");
       //  printf("------\n%c\n------",pPtr->parent->key);
      //last node now containing value info - Duplicate number of string/branch - First time 1.
      if(pPtr->children == NULL)
      {
         pPtr->children = TrieCreateNode('q', 1);
         pPtr->children->parent = pPtr;
      }
      else
      {
         //TrieValueNode_t *vPtr;
         pPtr->children->count = TrieCreateValueNode(data);
         pPtr->children->count->value ++;

         //if(pPtr->children->count->value > 2)
         //printf("Duplicates: %s -> %d\n\n",key, pPtr->children->count->value);
         //printf("------\n%d\n------",pPtr->children->value);
         //printf("Found same path but next pointer isn't NULL\nChildren -> value = Number of occurences of path\n");
      }
      return;
   }

   while (*key != '\0')
   {
      if (*key == pTrav->key)
      {
         key++;
#ifdef DEBUG
         printf("Traversing child: [%c]\n", pTrav->children->key);
#endif
         pTrav = pTrav->children;
      }
      else
         break;
   }

   while (pTrav->next)
   {
      if (*key == pTrav->next->key)
      {
         key++;
         TrieAdd(&(pTrav->next), key, 0);
         return;
      }
      pTrav = pTrav->next;
   }

   if (*key)
   {
      pTrav->next = TrieCreateNode(*key, 0xffffffff);
   }
   else
   {
      pTrav->next = TrieCreateNode(*key, 0);
   }

   pTrav->next->parent = pTrav->parent;
   pTrav->next->prev = pTrav;

#ifdef DEBUG
   printf("Inserting [%c] as neighbour of [%c] \n", pTrav->next->key, pTrav->key);
#endif

   if (!(*key))
      return;

   key++;

   for (pTrav = pTrav->next; *key; pTrav = pTrav->children)
   {
      pTrav->children = TrieCreateNode(*key, 0xffffffff);
      pTrav->children->parent = pTrav;
#ifdef DEBUG
      printf("Inserting: [%c]\n", pTrav->children->key);
#endif
      key++;
   }

//last node
   pTrav->children = TrieCreateNode('\0', 0);
   pTrav->children->parent = pTrav;
#ifdef DEBUG
   printf("Inserting: [%c]\n", pTrav->children->key);
#endif
   return;
}

trieNode_t *TrieSearch(trieNode_t *root, const char *key)
{
   trieNode_t *level = root;
   trieNode_t *pPtr = NULL;

   int lvl = 0;
   while (1)
   {
      trieNode_t *found = NULL;
      trieNode_t *curr;

      for (curr = level; curr != NULL; curr = curr->next)
      {
         if (curr->key == *key)
         {
            found = curr;
            lvl++;
            break;
         }
      }

      if (found == NULL)
         return NULL;

      if (*key == '\0')
      {
         pPtr = curr;
         return pPtr;
      }

      level = found->children;
      key++;
   }
}

void TrieDestroy(trieNode_t *root)
{
   trieNode_t *tPtr = root;
   trieNode_t *tmp = root;

   while (tPtr)
   {
      while (tPtr->children)
         tPtr = tPtr->children;

      if (tPtr->prev && tPtr->next)
      {
         tmp = tPtr;
         tPtr->next->prev = tPtr->prev;
         tPtr->prev->next = tPtr->next;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
      }
      else if (tPtr->prev && !(tPtr->next))
      {
         tmp = tPtr;
         tPtr->prev->next = NULL;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
      }
      else if (!(tPtr->prev) && tPtr->next)
      {
         tmp = tPtr;
         tPtr->parent->children = tPtr->next;
         tPtr->next->prev = NULL;
         tPtr = tPtr->next;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
      }
      else
      {
         tmp = tPtr;
         if (tPtr->parent == NULL)
         {
            /*Root*/
            free(tmp);
            return;
         }
         tPtr = tPtr->parent;
         tPtr->children = NULL;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
      }
   }
}




/*    Trie Remove function
   not useful in this implementation
   TreeDestroy, frees memory        

void TrieRemove(trieNode_t **root, char *key)
{
   trieNode_t *tPtr = NULL;
   trieNode_t *tmp = NULL;

   if (NULL == *root || NULL == key)
      return;

   tPtr = TrieSearch((*root)->children, key);

   if (NULL == tPtr)
   {
      printf("Key [%s] not found in trie\n", key);
      return;
   }

#ifdef DEBUG
   printf("Deleting key [%s] from trie\n", key);
#endif

   while (1)
   {
      if (tPtr->prev && tPtr->next)
      {
         tmp = tPtr;
         tPtr->next->prev = tPtr->prev;
         tPtr->prev->next = tPtr->next;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
         break;
      }
      else if (tPtr->prev && !(tPtr->next))
      {
         tmp = tPtr;
         tPtr->prev->next = NULL;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
         break;
      }
      else if (!(tPtr->prev) && tPtr->next)
      {
         tmp = tPtr;
         tPtr->next->prev = NULL;
         tPtr->parent->children = tPtr->next;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
         break;
      }
      else
      {
         tmp = tPtr;
         tPtr = tPtr->parent;
         tPtr->children = NULL;
#ifdef DEBUG
         printf("Deleted [%c] \n", tmp->key);
#endif
         free(tmp);
      }
   }

#ifdef DEBUG
   printf("Deleted key [%s] from trie\n", key);
#endif
}
                                                */