/*triedriver.c*/
/*
 * To Compile : gcc -o trie trie.c triedriver.c
 * To run: ./trie
 */
#include <stdio.h>
#include <stdlib.h>
#include "trie.h"

int main()
{
    trieNode_t *root;
    printf("Trie Example\n");
    
    /*Create a trie*/
    TrieCreate(&root);
    
    TrieAdd(&root, "01110011010", 1);
    TrieAdd(&root, "01110011001", 2);
    TrieAdd(&root, "11110011011", 3);
    TrieAdd(&root, "01110111010", 5);
    TrieRemove(&root, "01110011010");
    TrieAdd(&root, "00000000000", 6);
    TrieRemove(&root, "01110011010");
    TrieAdd(&root, "01110010110", 6);

    /*Destroy the trie*/
    TrieDestroy(root);
}

