
/*			  KD Tree Implementation
 *			       by Steven Wu
 *				3-29-2011
 * 	  
 *	An ideal tree is calculated with a given # of levels.
 *	The tree is built iteratively, and all N points are used to fill the tree.	
 *	After the tree is built, a new point p0 is placed in the grid.
 *	A recursive function is used to find the nearest point.
 *
 */
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <cstdlib>
#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>

static int dimsize = 2000;

struct Point {
    double x;
    double y;
    struct Point *left, *right;
};

/********************************************************
*	Return int value 2^level			*
*	Used to calculate  number of N for ideal tree	*
*	and # times to run iterations at each level	*
*********************************************************/
int twotothe( int level ) { 		
// a function to return int value 2^level
    int N =1;
    while ( level > 0 ) {
	N *=2;
	level--;
    }
    return (N);
}


/********************************************************
*    Print specified portions of array of struct Point	*
*********************************************************/
void printArray( int i_start, int i_end, struct Point arr[] ) {
    int i;
    for( i = i_start; i < i_end; i++ ) {
	printf(" %5d   %6.3f    %6.3f\n", i, arr[i].x, arr[i].y );
    }
}


/********************************************************************************
* 		Place N random Point(s) within dimsize box 			*
*********************************************************************************/

struct Point randoms( int N, struct Point *arr ) { 
    int i; 				    
 //   srand(time(NULL));
    int seed;
    seed = 11.7777;
    srand(seed);

    for( i = 0; i < N; i++) {
	arr[i].x = (double)dimsize * ( rand() % RAND_MAX  ) / RAND_MAX;
	arr[i].y = (double)dimsize * ( rand() % RAND_MAX  ) / RAND_MAX;
    }
    return *arr;
}

/********************************************************************************
*  		"Bubble Sort" is an O( N^2 ) Sorting Routine.			*
*  		Sort by x- or y- direction as determined by "dir" 		*
*   										*
*   It is seen to be very slow at large N ( depth 18 is N = 524,287 ) 		*
*   Computation time could be greatly reduced by implementing CUDA 		*
*   for large N.  For smaller N, sorting on CPU may be quicker			*
*   by not transferring data back and forth between CPU and GPU.		*
*										*
*   Empirically, sorting for depth = 14 (N ~ 32000) is acceptable (~7sec).	*
*   It would be most efficient to use a CUDA routine for depth > 14.		*
*   for N = 65535 (depth=15), computation time is (~1 min 10 sec)		*
*********************************************************************************/
struct Point bubbleSort( int i_start, int i_end, struct Point *array, int dir ) {	
    double tempx, tempy;
    int i, j;

    if( dir % 2 == 0 ) { 					//X-sort
	for ( i = i_start; i < i_end; i++ ) {
   	    for ( j = i_start; j < i_end-i-1+i_start; j++ ) {
		if( array[j].x > array[j+1].x )  {
		    tempx = array[j].x;
		    tempy = array[j].y;
		    array[j].x = array[j+1].x;
		    array[j].y = array[j+1].y;
		    array[j+1].x = tempx;
		    array[j+1].y = tempy;
	    	}
	    }
    	}
    }
    else {
    	for ( i = i_start; i < i_end; i++ ) {			//Y-sort
   	    for ( j = i_start; j < i_end-i-1+i_start; j++ ) {
	    	if( array[j].y > array[j+1].y )  {
		    tempx = array[j].x;
		    tempy = array[j].y;
		    array[j].x = array[j+1].x;
		    array[j].y = array[j+1].y;
		    array[j+1].x = tempx;
		    array[j+1].y = tempy;
	    	}
	    }
    	}
    }
    return *array;
}


float fillArrays( int i_start, int i_end, struct Point *a, float *tempx, float *tempy ) {
    // transfer struct members into temp arrays & reset index
    for( int i = i_start; i < i_end; i++ ) {
	tempx[i] = (a+i)->x;
	tempy[i] = (a+i)->y;

//	std::cout << i << "  " << tempx[i] << "  " << tempy[i] << std::endl;
    }
    return *tempx, *tempy;
}

int resetIndex( int i_start, int i_end, int *index ) {
    for( int i = i_start; i < i_end; i++ ) {
	index[i] = i;
    }
    return *index;
}


struct Point updateStruct( int i_start, int i_end, struct Point *a, float *tempx, 
			   float *tempy, int *index, int level ) {
    //update struct with sorted values in x or y dir
    if( level % 2 == 0 ) {
	for(int i = i_start; i < i_end; i++)  {
	    (a+i)->x = tempx[ i ];
	    (a+i)->y = tempy[ index[i] ];
	}
    }
    else if( level % 2 == 1 ) {		
	for(int i = i_start; i < i_end; i++)  {
	    (a+i)->x = tempx[ index[i] ];
	    (a+i)->y = tempy[ i ];
	}
    }
    return *a;
}

/************************************************************************
*	BINARY TREE OPERATIONS						*
*	http://cslibrary.stanford.edu/110/BinaryTrees.html 		*
*									*
* 	Original code author: Nick Parlante.				*
*  	Modified by Steven Wu.						*
*************************************************************************/			

/*-----------------------------------------------------------------------
 Helper function that allocates a new node 
 with the given data and NULL left and right Pointers.
------------------------------------------------------------------------*/
struct Point* NewNode( struct Point *p ) {						
    struct Point* parent;
    parent = (struct Point*)malloc(sizeof(struct Point));

    parent = p;
    parent->left = NULL;
    parent->right = NULL;

    return(parent);
}

/*-----------------------------------------------------------------------
 Given a binary search tree, and x & y coords, the new node is placed 
 in the tree correctly by alternating x- or y- dir according to depth.
 Starts at root with x condition, is (iterated) using "tick."		
-------------------------------------------------------------------------*/

struct Point *insert( int tick, struct Point* parent, struct Point* child ) { 				

    if (parent == NULL )		//  1.  If the tree is empty,
    	return  NewNode(child) ;	//      return a new, single node
 
    else if( tick % 2 == 1 )  {		//  2.  Otherwise, recur down the tree

     	if ( child->x <= parent->x )	
	    parent->left = insert( tick + 1, parent->left, child );
    	else 	
	    parent->right = insert( tick + 1, parent->right, child );

    	return(parent); // return the (unchanged) node pointer
    }
    else if( tick % 2 == 0 )  {

     	if ( child->y <= parent->y )	
	    parent->left = insert( tick + 1, parent->left, child );
    	else 
	    parent->right = insert( tick + 1, parent->right, child );

    	return(parent); // return the (unchanged) node pointer
    }
    return;
}

	
/*****************************************************************************
*	A recursive function to find the leaf closest to the new point p.    * 
******************************************************************************/

struct Point simple_fn( int tick, struct Point p, struct Point *parent ) {

    if( parent->left == NULL && parent->right == NULL )    //return if leaf	
	return *parent;					

    else if( tick % 2 == 0)  {			// x- dim comparison
	tick++;
	if( p.x <= parent->x )
	    parent = parent->left;
    	else 
	    parent = parent->right;
	return simple_fn( tick, p, parent );	//call function, with tick+1
    }

    else if( tick % 2 == 1 )  {			// y- dim comparison
	tick++;
	if( p.y <= parent->y )
	    parent = parent->left;
    	else 
	    parent = parent->right;
	return simple_fn( tick, p, parent );
    }
    return;
}



  //======================================================================//
 //	main() starts here					 	 //
//======================================================================//

main()  
{
    int i;				// array index
    int N = 1; 				// Number of N for the level = 0 case
    int level_tot = 26;			// Enter the total number of levels
    int level = level_tot;              // Initialize level counter
    int i_start, i_end;			// Index start/end of sub-arrays
    int median;				// Median value chosen multiple times 
    int tick;				// Used to ensure leaf runs through 
					// conditions starting at root
    double distx, disty, distance;

   while ( level > 0 )  {		// For a given level,
	N += twotothe( level );         // calculate N needed for
	level--;			// balanced tree
    }

    int N_tot = N;			// Retain N value for after tree is built

    struct Point *a;  			// one Point of x and y coords
    a = (struct Point*)malloc( N * sizeof( *a ) );     // Allocate memory for N Points

    struct Point* root;
    root = NULL;

    float *tempx, *tempy;
    tempx = (float*)malloc( N * sizeof( *tempx) );
    tempy = (float*)malloc( N * sizeof( *tempy) );

    int *index;
    index = (int*)malloc( N * sizeof( *index ) ); 


    randoms( N, a );    		// place N random pts in box

//  printArray(0, N_tot, a);

    level = 0;				// Start at root


/********************************************************************************
*   Sort initial array in x- dir, find median, and place coords in root.	*
*   Split into two arrays at the median, 					*
*   sort each array in y- dir and place each median in tree.			*
*   Switching off between x- and y- conditions,					*
*   repeat until while() condition is met, and place coords correctly in tree.  *
*********************************************************************************/
	
    while( N >=3 ) {

//	Uncomment print statements to check points at every level
	printf("\nLEVEL %d\n", level);

   	for( i = 0; i < twotothe(level); i++ )  {         

	    i_start = i * ( N + 1 );
	    i_end   = i_start + N;    
	    printf("i = %d to  i < %d:\n", i_start, i_end);

	    fillArrays( i_start, i_end, a, tempx, tempy );
	    resetIndex( i_start, i_end, index );

//	    bubbleSort( i_start, i_end, a, level );
 
	    if( level % 2 == 0 )	
		thrust::stable_sort_by_key( tempx + i_start, tempx + i_end, index + i_start); 
	    else 	
		thrust::stable_sort_by_key( tempy + i_start, tempy + i_end, index + i_start); 

    	    updateStruct( i_start, i_end, a, tempx, tempy, index, level );
	    printArray( i_start, i_end, a );			

	    median = N/2+i*(N+1);
	    printf("median = %d\n", median );

	    tick = 1;		// reset the ticker before each call to insert()
	    root = insert( tick, root, ( a + median ) );
	}

    	level++;
    	N/=2;
    }


/************************************************************************
*		Find the closest point to p0.				*
*	If the distance to a leaf's parent is closer than the leaf,	*
*	this function call does not return the closest point		*
*************************************************************************/


    struct Point p0;

    p0.x = 3;		// Arbitrarily chosen points within dimsize.
    p0.y = 18;		// Round numbers avoid running into same value in grid.

    for( i = 0; i < N_tot; i++ ) {
	distx = p0.x - (a+i)->x;
	distx *= distx;
	disty = p0.y - (a+i)->y;
	disty *= disty;
	distance = sqrt(distx + disty);
	printf("%d  ( %f   %f )     %f\n", i, (a+i)->x, (a+i)->y, distance);
    }
    printf("**************************************\n");    
	
    level = 0;

    *root = simple_fn( level, p0, root ) ;
    printf("Closest point to p0: %f %f\n", root->x, root->y );

    printf("N = %d\n", N_tot);       	 
    printf("total levels in tree: %d\n", level_tot);

    free( a );
    free( tempx );
    free( tempy );
    free( index );

    return 0;

}



