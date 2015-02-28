#pragma once

#include <vector>
#include <stdexcept>
#include <math.h>
#include <cstring>

#define kMaxPointsPerCell 9

class RangeList;
class ScallopedRegion;

double timeInSeconds(void);

class RNG {
private:
	/* Period parameters */
	static const long N = 624;
	static const long M = 397;
	static const unsigned long MATRIX_A = 0x9908b0dfUL;   /* constant vector a */
	static const unsigned long UPPER_MASK = 0x80000000UL; /* most significant w-r bits */
	static const unsigned long LOWER_MASK = 0x7fffffffUL; /* least significant r bits */

private:
	unsigned long mt[N]; /* the array for the state vector  */
	int mti;

public:
	RNG(unsigned long seed=5489UL);
	RNG(unsigned long *init_key, int key_length);

	void seed(unsigned long seed);

		/* generates a random number on [0,0xffffffff]-interval */
	unsigned long getInt32();
		/* generates a random number on [0,0x7fffffff]-interval */
	long getInt31();
		/* generates a random number on [0,1]-real-interval */
	double getDoubleLR();
	float getFloatLR();
		/* generates a random number on [0,1)-real-interval */
	double getDoubleL();
	float getFloatL();
		/* generates a random number on (0,1)-real-interval */
	double getDouble();
	float getFloat();
};

class Vec2 {
public:
	Vec2() {};
	Vec2(float _x, float _y) : x(_x), y(_y) {};

	float x,y;

	float length() { return sqrt(x*x + y*y); }

	bool operator ==(const Vec2 &b) const { return x==b.x && y==b.y; }
	Vec2 operator +(Vec2 b) { return Vec2(x+b.x, y+b.y); }
	Vec2 operator -(Vec2 b) { return Vec2(x-b.x, y-b.y); }
	Vec2 operator *(Vec2 b) { return Vec2(x*b.x, y*b.y); }
	Vec2 operator /(Vec2 b) { return Vec2(x/b.x, y*b.y); }

	Vec2 operator +(float n) { return Vec2(x+n, y+n); }
	Vec2 operator -(float n) { return Vec2(x-n, y-n); }
	Vec2 operator *(float n) { return Vec2(x*n, y*n); }
	Vec2 operator /(float n) { return Vec2(x/n, y*n); }

	Vec2 &operator +=(Vec2 b) { x+=b.x; y+=b.y; return *this; }
	Vec2 &operator -=(Vec2 b) { x-=b.x; y-=b.y; return *this; }
	Vec2 &operator *=(Vec2 b) { x*=b.x; y*=b.y; return *this; }
	Vec2 &operator /=(Vec2 b) { x/=b.x; y/=b.y; return *this; }

	Vec2 &operator +=(float n) { x+=n; y+=n; return *this; }
	Vec2 &operator -=(float n) { x-=n; y-=n; return *this; }
	Vec2 &operator *=(float n) { x*=n; y*=n; return *this; }
	Vec2 &operator /=(float n) { x/=n; y/=n; return *this; }
};

class PDSampler {
protected:
	RNG m_rng;
	std::vector<int> m_neighbors;

	int (*m_grid)[kMaxPointsPerCell];
	int m_gridSize;
	float m_gridCellSize;

public:
	std::vector<Vec2> points;
	float radius;
	bool isTiled;

public:
	PDSampler(float radius, bool isTiled, bool usesGrid=true);
	virtual ~PDSampler() { };

	//

	bool pointInDomain(Vec2 &a);

		// return shortest distance between _a_
		// and _b_ (accounting for tiling)
	float getDistanceSquared(Vec2 &a, Vec2 &b) { Vec2 v = getTiled(b-a); return v.x*v.x + v.y*v.y; }
	float getDistance(Vec2 &a, Vec2 &b) { return sqrt(getDistanceSquared(a, b)); }

		// generate a random point in square
	Vec2 randomPoint();

		// return tiled coordinates of _v_
	Vec2 getTiled(Vec2 v);

		// return grid x,y for point
	void getGridXY(Vec2 &v, int *gx_out, int *gy_out);

		// add _pt_ to point list and grid
	void addPoint(Vec2 pt);

		// populate m_neighbors with list of
		// all points within _radius_ of _pt_
		// and return number of such points
	int findNeighbors(Vec2 &pt, float radius);

		// return distance to closest neighbor within _radius_
	float findClosestNeighbor(Vec2 &pt, float radius);

		// find available angle ranges on boundary for candidate
		// by subtracting occluded neighbor ranges from _rl_
	void findNeighborRanges(int index, RangeList &rl);

		// extend point set by boundary sampling until domain is
		// full
	void maximize();

		// apply one step of Lloyd relaxation
	void relax();

	//

	virtual void complete() = 0;
};

///

class DartThrowing : public PDSampler {
private:
	int m_minMaxThrows, m_maxThrowsMult;

public:
	DartThrowing(float radius, bool isTiled, int minMaxThrows, int maxThrowsMult);

	virtual void complete();
};


class BestCandidate : public PDSampler {
private:
	int m_multiplier, m_N;

public:
	BestCandidate(float radius, bool isTiled, int multiplier);

	virtual void complete();
};

///

class BoundarySampler : public PDSampler {
public:
	BoundarySampler(float radius, bool isTiled) : PDSampler(radius, isTiled) {};

	virtual void complete();
};

///

class PureSampler : public PDSampler {
public:
	PureSampler(float radius) : PDSampler(radius, true) {};

	virtual void complete();
};

///

class LinearPureSampler : public PDSampler {
public:
	LinearPureSampler(float radius) : PDSampler(radius, true) {};

	virtual void complete();
};

///

class PenroseSampler : public PDSampler {
public:
	PenroseSampler(float radius) : PDSampler(radius, false, false) {};

	virtual void complete();
};

///

class UniformSampler : public PDSampler {
public:
	UniformSampler(float radius) : PDSampler(radius, false, false) {};

	virtual void complete();
};

template <class T>
class WDPDF_Node
{
private:
	bool m_mark;

public:
	WDPDF_Node<T> *parent, *left, *right;
	T key;
	float weight, sumWeights;

public:
	WDPDF_Node(T key_, float weight_, WDPDF_Node<T> *parent_);
	~WDPDF_Node();

	WDPDF_Node<T> *sibling() { return this==parent->left?parent->right:parent->left; }

	void markRed() { m_mark = true; }
	void markBlack() { m_mark = false; }
	bool isRed() { return m_mark; }
	bool isBlack() { return !m_mark; }
	bool leftIsBlack() { return !left || left->isBlack(); }
	bool rightIsBlack() { return !right || right->isBlack(); }
	bool leftIsRed() { return !leftIsBlack(); }
	bool rightIsRed() { return !rightIsBlack(); }
	void setSum() { sumWeights = weight + (left?left->sumWeights:0) + (right?right->sumWeights:0); }
};

template <class T>
class WeightedDiscretePDF
{
private:
	WDPDF_Node<T> *m_root;

public:
	WeightedDiscretePDF();
	~WeightedDiscretePDF();

	void insert(T item, float weight);
	void update(T item, float newWeight);
	void remove(T item);
	bool inTree(T item);

		/* pick a tree element according to its
		 * weight. p should be in [0,1).
		 */
	T choose(float p);

private:
	WDPDF_Node<T> **lookup(T item, WDPDF_Node<T> **parent_out);
	void split(WDPDF_Node<T> *node);
	void rotate(WDPDF_Node<T> *node);
	void lengthen(WDPDF_Node<T> *node);
	void propogateSumsUp(WDPDF_Node<T> *n);
};


template <class T>
WDPDF_Node<T>::WDPDF_Node(T key_, float weight_, WDPDF_Node<T> *parent_)
{
	m_mark = false;

	key = key_;
	weight = weight_;
	sumWeights = 0;
	left = right = 0;
	parent = parent_;
}

template <class T>
WDPDF_Node<T>::~WDPDF_Node()
{
	if (left) delete left;
	if (right) delete right;
}

//

template <class T>
WeightedDiscretePDF<T>::WeightedDiscretePDF()
{
	m_root = 0;
}

template <class T>
WeightedDiscretePDF<T>::~WeightedDiscretePDF()
{
	if (m_root) delete m_root;
}

template <class T>
void WeightedDiscretePDF<T>::insert(T item, float weight)
{
	WDPDF_Node<T> *p=0, *n=m_root;

	while (n) {
		if (n->leftIsRed() && n->rightIsRed())
			split(n);

		p = n;
		if (n->key==item) {
			throw std::domain_error("insert: argument(item) already in tree");
		} else {
			n = (item<n->key)?n->left:n->right;
		}
	}

	n = new WDPDF_Node<T>(item, weight, p);

	if (!p) {
		m_root = n;
	} else {
		if (item<p->key) {
			p->left = n;
		} else {
			p->right = n;
		}

		split(n);
	}

	propogateSumsUp(n);
}

template <class T>
void WeightedDiscretePDF<T>::remove(T item)
{
	WDPDF_Node<T> **np = lookup(item, 0);
	WDPDF_Node<T> *child, *n = *np;

	if (!n) {
		throw std::domain_error("remove: argument(item) not in tree");
	} else {
		if (n->left) {
			WDPDF_Node<T> **leftMaxp = &n->left;

			while ((*leftMaxp)->right)
				leftMaxp = &(*leftMaxp)->right;

			n->key = (*leftMaxp)->key;
			n->weight = (*leftMaxp)->weight;

			np = leftMaxp;
			n = *np;
		}

			// node now has at most one child

		child = n->left?n->left:n->right;
		*np = child;

		if (child) {
			child->parent = n->parent;

			if (n->isBlack()) {
				lengthen(child);
			}
		}

		propogateSumsUp(n->parent);

		n->left = n->right = 0;
		delete n;
	}
}

template <class T>
void WeightedDiscretePDF<T>::update(T item, float weight)
{
	WDPDF_Node<T> *n = *lookup(item, 0);

	if (!n) {
		throw std::domain_error("update: argument(item) not in tree");
	} else {
		float delta = weight - n->weight;
		n->weight = weight;

		for (; n; n=n->parent) {
			n->sumWeights += delta;
		}
	}
}

template <class T>
T WeightedDiscretePDF<T>::choose(float p)
{
	if (p<0.0 || p>=1.0) {
		throw std::domain_error("choose: argument(p) outside valid range");
	} else if (!m_root) {
		throw std::logic_error("choose: choose() called on empty tree");
	} else {
		float w = m_root->sumWeights * p;
		WDPDF_Node<T> *n = m_root;

		while (1) {
			if (n->left) {
				if (w<n->left->sumWeights) {
					n = n->left;
					continue;
				} else {
					w -= n->left->sumWeights;
				}
			}
			if (w<n->weight || !n->right) {
				break; // !n->right condition shouldn't be necessary, just sanity check
			}
			w -= n->weight;
			n = n->right;
		}

		return n->key;
	}
}

template <class T>
bool WeightedDiscretePDF<T>::inTree(T item)
{
	WDPDF_Node<T> *n = *lookup(item, 0);

	return !!n;
}

//

template <class T>
WDPDF_Node<T> **WeightedDiscretePDF<T>::lookup(T item, WDPDF_Node<T> **parent_out)
{
	WDPDF_Node<T> *n, *p=0, **np=&m_root;

	while ((n = *np)) {
		if (n->key==item) {
			break;
		} else {
			p = n;
			if (item<n->key) {
				np = &n->left;
			} else {
				np = &n->right;
			}
		}
	}

	if (parent_out)
		*parent_out = p;
	return np;
}

template <class T>
void WeightedDiscretePDF<T>::split(WDPDF_Node<T> *n)
{
	if (n->left) n->left->markBlack();
	if (n->right) n->right->markBlack();

	if (n->parent) {
		WDPDF_Node<T> *p = n->parent;

		n->markRed();

		if (p->isRed()) {
			p->parent->markRed();

				// not same direction
			if (!(	(n==p->left && p==p->parent->left) ||
					(n==p->right && p==p->parent->right))) {
			  rotate(n);
			  p = n;
			}

			rotate(p);
			p->markBlack();
		}
	}
}

template <class T>
void WeightedDiscretePDF<T>::rotate(WDPDF_Node<T> *n)
{
	WDPDF_Node<T> *p=n->parent, *pp=p->parent;

	n->parent = pp;
	p->parent = n;

	if (n==p->left) {
		p->left = n->right;
		n->right = p;
		if (p->left) p->left->parent = p;
	} else {
		p->right = n->left;
		n->left = p;
		if (p->right) p->right->parent = p;
	}

	n->setSum();
	p->setSum();

	if (!pp) {
		m_root = n;
	} else {
		if (p==pp->left) {
			pp->left = n;
		} else {
			pp->right = n;
		}
	}
}

template <class T>
void WeightedDiscretePDF<T>::lengthen(WDPDF_Node<T> *n)
{
	if (n->isRed()) {
		n->markBlack();
	} else if (n->parent) {
		WDPDF_Node<T> *sibling = n->sibling();

		if (sibling && sibling->isRed()) {
			n->parent->markRed();
			sibling->markBlack();

			rotate(sibling); // node sibling is now old sibling child, must be black
			sibling = n->sibling();
		}

		// sibling is black

		if (!sibling) {
			lengthen(n->parent);
		} else if (sibling->leftIsBlack() && sibling->rightIsBlack()) {
			if (n->parent->isBlack()) {
				sibling->markRed();
				lengthen(n->parent);
			} else {
				sibling->markRed();
				n->parent->markBlack();
			}
		} else {
			if (n==n->parent->left && sibling->rightIsBlack()) {
				rotate(sibling->left); // sibling->left must be red
				sibling->markRed();
				sibling->parent->markBlack();
				sibling = sibling->parent;
			} else if (n==n->parent->right && sibling->leftIsBlack()) {
				rotate(sibling->right); // sibling->right must be red
				sibling->markRed();
				sibling->parent->markBlack();
				sibling = sibling->parent;
			}

			// sibling is black, and sibling's far child is red

			rotate(sibling);
			if (n->parent->isRed()) sibling->markRed();
			sibling->left->markBlack();
			sibling->right->markBlack();
		}
	}
}

template <class T>
void WeightedDiscretePDF<T>::propogateSumsUp(WDPDF_Node<T> *n)
{
	for (; n; n=n->parent)
		n->setSum();
}

typedef struct _RangeEntry {
	float min, max;
} RangeEntry;

class RangeList {
public:
	RangeEntry *ranges;
	int numRanges, rangesSize;

public:
	RangeList(float min, float max);
	~RangeList();

	void reset(float min, float max);

	void print();

	void subtract(float min, float max);

private:
	void deleteRange(int pos);
	void insertRange(int pos, float min, float max);
};

typedef struct {
	Vec2 P;
	float r, sign, d, theta, integralAtStart;
	float rSqrd, dSqrd;
} ArcData;

class ScallopedSector
{
public:
	Vec2 P;
	float a1, a2, area;

	ArcData arcs[2];

public:
	ScallopedSector(Vec2 &_Pt, float _a1, float _a2, Vec2 &P1, float r1, float sign1, Vec2 &P2, float r2, float sign2);

	float calcAreaToAngle(float angle);
	float calcAngleForArea(float area, RNG &rng);
	Vec2 sample(RNG &rng);

	float distToCurve(float angle, int index);

	void subtractDisk(Vec2 &C, float r, std::vector<ScallopedSector> *regions);

private:
	float canonizeAngle(float angle);

	void distToCircle(float angle, Vec2 &C, float r, float *d1_out, float *d2_out);
};

class ScallopedRegion
{
public:
	std::vector<ScallopedSector> *regions;
	float minArea;
	float area;

public:
	ScallopedRegion(Vec2 &P, float r1, float r2, float minArea=.00000001);
	~ScallopedRegion();

	bool isEmpty() { return regions->size()==0; }
	void subtractDisk(Vec2 C, float r);

	Vec2 sample(RNG &rng);
};


#ifndef QUASISAMPLER_PROTOTYPE_H
#define QUASISAMPLER_PROTOTYPE_H

#include <math.h>
#include <vector>

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define LUT_SIZE 21 // Number of Importance Index entries in the Lookup table.
#define NUM_STRUCT_INDEX_BITS 6 // Number of significant bits taken from F-Code.

#define GOLDEN_RATIO PHI // Phi is the Golden Ratio.
#define PHI     1.6180339887498948482045868343656 // ( 1 + sqrt(5) ) / 2
#define PHI2    2.6180339887498948482045868343656 // Phi squared
#define LOG_PHI 0.48121182505960347 // log(Phi)
#define SQRT5   2.2360679774997896964091736687313 // sqrt(5.0)

// Two-bit sequences.
#define B00 0
#define B10 1
#define B01 2

/// The six tile types.
enum TileType {
  TileTypeA,TileTypeB,TileTypeC,TileTypeD,TileTypeE,TileTypeF
};

/// Simple 2D point and vector type.
class Point2D
{
public:
  double x,y;

  Point2D(){};
  Point2D(const double x, const double y) { this->x=x; this->y=y; }
  Point2D(const double vect[2]) { x=vect[0]; y=vect[1]; }

  Point2D operator+(const Point2D& pt) const{ return Point2D(x+pt.x,y+pt.y); }
  Point2D operator-(const Point2D& pt) const{ return Point2D(x-pt.x,y-pt.y); }
  Point2D operator*(double factor) const{ return Point2D(x*factor,y*factor); }
  Point2D operator/(double factor) const{ return Point2D(x/factor,y/factor); }

  /// Returns the squared distance to the origin, or the squared length of a vector.
  double d2() const { return x*x+y*y; }
};

/// This is a base class that implements the Quasi-Sampler importance sampling
/// system, as presented in the paper :
/// "Fast Hierarchical Importance Sampling with Blue Noise Properties",
/// by Victor Ostromoukhov, Charles Donohue and Pierre-Marc Jodoin,
/// to be presented at SIGGRAPH 2004.
/// This is a pure-virtual class, and you must implement the "getImportanceAt()" function
/// in order to use the sampling system.
/// The mechanics of the system can be observed in the given source code.
class Quasisampler
{

protected:

  //
  // Static tables.
  //


  /// Fibonacci sequence (first 32 numbers).
  static const unsigned fiboTable[32]; // defined at end of file.

  /// Unit vectors rotated around origin, in \f$ \frac{\pi}{10} \f$ increments,
  /// counter-clockwise. 0 = North.
  /// This table can be used to accelerate the trigonomic operations within the tile
  /// subdivision process, since all angles can only take these values.
  static const Point2D vvect[20]; // defined at end of file.

  /// Pre-calculated correction vectors lookup table.
  /// These are available in ASCII format on the web-site.
  static const double lut[LUT_SIZE][21][2]; // defined at end of file.



  //
  //  Static functions.
  //

  /// Fibonacci number at a given position.
  /// The value returned is \f$ F_i = F_{i-1} + F_{i-2}  \f$.
  static unsigned fibonacci(unsigned i)
  {
    if (i<1) return 1;
    if (i<=32) return fiboTable[i-1]; // pre-calculated.
    return fibonacci(i-1)+fibonacci(i-2);
  }

  /// Returns the required level of subdivision for a given importance value.
  /// The value returned is \f$ \lceil{\log_{\phi^2}(importance)}\rceil \f$,
  /// where \f$ \phi=\frac{1 + {\sqrt{5}}}{2}\f$  is the Golden Ratio.
  static unsigned getReqSubdivisionLevel( unsigned importance )
  {
    if (importance==0) return 0;
    unsigned nbits = (unsigned)(log( (double)importance*SQRT5 + 1.0 ) / LOG_PHI) - 1;
    if (nbits<1) nbits = 1;
    return (unsigned)ceil(0.5*nbits);
  }

  /// Returns the decimal value of an F-Code, over a given number of bits.
  /// The value returned is \f$ \sum_{j=2}^{m} b_{j} F_{j} \f$.
  static unsigned calcFCodeValue(unsigned bitsequence,unsigned nbits)
  {
    unsigned i_s = 0;
    for (unsigned i=0; i<nbits ;i++ )
    {
      if ( bitsequence & ( 1u<<(nbits-i-1) ) ) i_s += fibonacci(i+2);
    }
    return i_s;
  }


  /// Returns the Structural Index (i_s) for a given F-Code.
  static unsigned calcStructuralIndex(unsigned bitsequence)
  {
    return calcFCodeValue(bitsequence,NUM_STRUCT_INDEX_BITS);
  }

  /// Returns the Importance Index (i_v) for a given importance value.
  /// The value returned is \f$ \lfloor n \cdot ({\log_{\phi^2} \sqrt{5} \cdot x}) ~ {\bf mod} ~ 1 \rfloor \f$.
  static unsigned calcImportanceIndex( unsigned importance )
  {
    double t = log(1.0 + sqrt(5.0)*importance) / log(PHI2);
    t -= floor(t); // modulo 1.0
    return (unsigned)(LUT_SIZE*t);
  }


  /// Fetches the appropriate vector from the lookup table.
  static Point2D calcDisplacementVector(unsigned importance, unsigned f_code, int dir)
  {
    unsigned i_s = calcStructuralIndex(f_code);
    unsigned i_v = calcImportanceIndex(importance);

    return
      vvect[dir]        * lut[i_v][i_s][0] + // u component
      vvect[(dir+5)%20] * lut[i_v][i_s][1] ; // v component
  }


  //
  // Inner classes.
  //


  /// Individual tile elements, which also serve as nodes for the tile subdivision tree.
  class TileNode
  {

    unsigned level; // Depth in the tree.
    int tileType; // Types A through F.
    int dir; // Tile orientation, 0=North, in Pi/10 increments, CCW.
    double scale;
    Point2D p1,p2,p3; // Three points of the triangle. Counter-clockwise.

    /// The F-Code binary sequence.
    unsigned f_code;

    // tiling tree structure
    TileNode* parent;
    unsigned parent_slot; // position in parent's list (needed for iterators)
    bool terminal; // true for leaf nodes
    std::vector<TileNode*> children;

  public:

    /// Builds a tile according to the given specifications.
    TileNode(
      TileNode* parent = NULL,
      int tileType = TileTypeF,
      Point2D refPt = Point2D(0,0),
      int dir = 15, // 15 = East.
      unsigned newbits = 0,
      int parent_slot = 0,
      double scale = 1.0)
    {
      this->parent = parent;
      this->tileType = tileType;
      this->p1 = refPt;
      this->dir = dir%20;
      this->parent_slot = parent_slot;
      this->scale = scale;
      this->level = parent ? parent->level + 1 : 0; // Increment the level.


      // Build triangle, according to type.
      switch(tileType)
      {
      case TileTypeC:
      case TileTypeD:
        // "Skinny" triangles
        p2 = p1 + vvect[dir%20]*scale;
        p3 = p1 + vvect[(dir+4)%20]*(PHI*scale);
        break;
      case TileTypeE:
      case TileTypeF:
        // "Fat" triangles
        p2 = p1 + vvect[dir%20]*(PHI2*scale);
        p3 = p1 + vvect[(dir+2)%20]*(PHI*scale);
        break;
      default:
        // Pentagonal tiles (triangle undefined)
        p2 = p1 + vvect[dir%20]*scale;
        p3 = p1 + vvect[(dir+5)%20]*scale;
      }

      // Append 2 new bits to the F-Code.
      if (parent)
        f_code = (parent->f_code<<2)^newbits;
      else
        f_code = newbits;

      // Set as leaf node
      terminal = true;
      children.clear();
    }

    /// Helper constructor.
    /// Creates an initial tile that is certain to contain the ROI.
    /// The starting tile is of type F (arbitrary).
    TileNode( double roi_width,  double roi_height)
    {
      double side = MAX(roi_width,roi_height);
      double scale = 2.0 * side;
      Point2D offset(PHI*PHI/2.0-0.25,0.125);
      *this = TileNode(NULL, TileTypeF,offset*-scale,15,0,0,scale);
    }

    ~TileNode() { collapse(); }

    /// Splits a tile according to the given subdivision rules.
    /// Please refer to the code for further details.
    void refine()
    {
      if (!terminal) return; // Can only subdivide leaf nodes.

      terminal=false; // The tile now has children.
      double newscale = scale / GOLDEN_RATIO; // The scale factor between levels is constant.

      switch(tileType)
      {

      // Each new tile is created using the following information:
      // A pointer to its parent, the type of the new tile (a through f),
      // the origin of the new tile, the change in orientation of the new tile with
      // respect to the parent's orientation, the two bits to be pre-pended to the F-Code,
      // the parent's slot (for traversal purposes), and the new linear scale of the tile,
      // which is always the parent's scale divided by the golden ratio.

      case TileTypeA:
        children.push_back(
          new TileNode(this, TileTypeB, p1, dir+0, B00, 0, newscale));
        break;

      case TileTypeB:
        children.push_back(
          new TileNode(this, TileTypeA, p1, dir+10, B00, 0, newscale) );
        break;

      case TileTypeC:
        children.push_back(
          new TileNode(this, TileTypeF, p3, dir+14, B00, 0, newscale) );
        children.push_back(
          new TileNode(this, TileTypeC, p2, dir+6, B10, 1, newscale) );
        children.push_back(
          new TileNode(this, TileTypeA, children[0]->p3, dir+1, B10, 2, newscale) );
        break;

      case TileTypeD:
        children.push_back(
          new TileNode(this, TileTypeE, p2, dir+6, B00, 0, newscale) );
        children.push_back(
          new TileNode(this, TileTypeD, children[0]->p3, dir+14, B10, 1, newscale) );
        break;

      case TileTypeE:
        children.push_back(
          new TileNode(this, TileTypeC, p3, dir+12, B10, 0, newscale) );
        children.push_back(
          new TileNode(this, TileTypeE, p2, dir+8 , B01, 1, newscale) );
        children.push_back(
          new TileNode(this, TileTypeF, p1, dir+0 , B00, 2, newscale) );
        children.push_back(
          new TileNode(this, TileTypeA, children[0]->p2, dir+7, B10, 3, newscale) );
        break;

      case TileTypeF:
        children.push_back(
          new TileNode(this, TileTypeF, p3, dir+12, B01, 0, newscale) );
        children.push_back(
          new TileNode(this, TileTypeE, children[0]->p3, dir+0, B00, 1, newscale) );
        children.push_back(
          new TileNode(this, TileTypeD, children[1]->p3, dir+8, B10, 2, newscale) );
        children.push_back(
          new TileNode(this, TileTypeA, children[0]->p3, dir+15, B01, 3, newscale) );
        break;
      }
    }

    /// Prunes the subdivision tree at this node.
    void collapse()
    {
      // Recursively prune the tree.
      for (unsigned i=0; i<children.size(); i++) delete children[i];
      terminal = true;
      children.clear();
    }

    /// Returns the next node of the tree, in depth-first traversal.
    /// Returns NULL if it is at the last node.
    TileNode* nextNode()
    {
      if (!terminal) return children[0];

      if (level == 0) return NULL; // single node case.

      if ( parent_slot < parent->children.size()-1 )
        return parent->children[parent_slot+1];

      // last child case
      TileNode* tmp = this;
      do
      {
        tmp = tmp->parent;
      }
      while ( (tmp->level != 0) && (tmp->parent_slot ==  tmp->parent->children.size()-1) );

      if (tmp->level == 0) return NULL; // last node
      return tmp->parent->children[tmp->parent_slot+1];

    }

    /// Returns the next closest leaf to a node.
    /// Returns NULL if it's the last leaf.
    TileNode* nextLeaf()
    {
      TileNode* tmp = this;
      do
      {
        tmp = tmp->nextNode();
        if ( !tmp ) return NULL;
        if ( tmp->terminal ) return tmp;
      }
      while (1);
    }

    // Public accessors

    Point2D getP1() const { return p1; }
    Point2D getP2() const { return p2; }
    Point2D getP3() const { return p3; }
    Point2D getCenter() const { return (p1+p2+p3)/3.0; }
    unsigned getFCode() const { return f_code; }
    bool isSamplingType() const {
      return ( (tileType == TileTypeA) || (tileType == TileTypeB) ); }
    unsigned getLevel() { return level; }
    bool isTerminal() const { return terminal; }
    TileNode* getParent() { return parent; }
    TileNode* getChild(unsigned i) {  return children[i]; }

    /// Obtains the correction vector from the lookup table,
    /// then scales and adds it to the reference point.
    Point2D getDisplacedSamplingPoint(unsigned importance)
    {
      return p1 + calcDisplacementVector(importance,f_code,dir) * scale;
    }

  }; // end of class TileNode.

  /// Leaf iterator for the tile subdivision tree.
  /// The traversal is made in a depth-first manner.
  /// Warning: This does not behave like STL style iterators.
  class TileLeafIterator
  {
    TileNode* shape;
  public:
    TileLeafIterator() { shape=NULL; }
    TileLeafIterator(TileNode* s ) { begin(s); }

    TileNode* operator*() { return shape; }
    TileNode* operator->() { return shape; }

    void begin(TileNode* s)
    {
      TileNode* tmp = s;
      while ( ! tmp->isTerminal() ) tmp = tmp->getChild(0); // find first leaf
      shape = tmp;
    }

    /// Subdivides the tile and moves to its 1st child.
    void refine()
    {
      shape->refine();
      shape = shape->getChild(0);
    }

    /// Prunes the subdivision tree.
    void collapse()
    {
      if (shape->getParent())
      {
        shape = shape->getParent();
        shape->collapse();
      }
    }

    /// Moves to the next node in the subdivision tree, in depth-first traversal.
    /// Returns false iff there is no such node.
    bool next()
    {
      TileNode* s = shape->nextLeaf();
      if (s)
      {
        shape = s;
        return true;
      }
      else
      {
        shape = s;
        return false;
      }
    }

    /// Checks if there is a next tile, in depth-first traversal.
    bool hasNext()
    {
      TileNode* s = shape->nextLeaf();
      if (s) return true;
      else return false;
    }
  };


  //
  // Instance members.
  //

  /// Root node of the tile subdivision tree.
  TileNode *root;

  /// Extents of the region of interest.
  double width, height;

  /// Protected constructor, which initializes the Region of Interest.
  Quasisampler(double width=0.0, double height=0.0)
  { this->width=width; this->height=height; root=NULL; }

  virtual ~Quasisampler() { if (root) delete root; }



  /// This is a helper function which constrains the incoming points
  /// to the region of interest.
  unsigned getImportanceAt_bounded(Point2D pt)
  {
    if (pt.x>=0 && pt.x<width && pt.y>=0 && pt.y<height)
      return getImportanceAt(pt);
    else
      return 0;
  }

  /// Subdivides all tiles down a level, a given number of times.
  void subdivideAll(int times=1)
  {
    if (!root) return;
    TileNode *tmp;
    for (int i=0;i<times;i++)
    {
      TileLeafIterator it(root);
      do {
        tmp = *it;
        it.next();
        tmp->refine();
      }
      while (*it);
    }
  }

  /// Generates the hierarchical structure.
  void buildAdaptiveSubdivision( unsigned minSubdivisionLevel = 6 )
  {
    root = new TileNode(width,height);

    // Since we are approximating the MAX within each tile by the values at
    // a few key points, we must provide a sufficiently dense initial
    // tiling. This would not be necessary with a more thorough scan of each
    // tile.
    subdivideAll(minSubdivisionLevel);

    TileLeafIterator it(root);
    TileNode *tmp;

    // Recursively subdivide all triangles until each triangle's
    // required level is reached.
    unsigned level;
    do {
      level = it->getLevel();

      if ( it->isSamplingType() ) // Sampling tiles are infinitesimal
      {
        if ( level < getReqSubdivisionLevel(getImportanceAt_bounded(it->getP1())) )
        {
          tmp = *it;
          tmp->refine();
        }
      }
      else
      {
        if (
          ( level < getReqSubdivisionLevel(getImportanceAt_bounded(it->getP1())) ) ||
          ( level < getReqSubdivisionLevel(getImportanceAt_bounded(it->getP2())) ) ||
          ( level < getReqSubdivisionLevel(getImportanceAt_bounded(it->getP3())) ) ||
          ( level < getReqSubdivisionLevel(getImportanceAt_bounded(it->getCenter())) )
          )
        {
          tmp = *it;
          tmp->refine();
        }
      }
    } while ( it.next() );
  }


  /// Collect the resulting point set.
  void collectPoints(
    std::vector<Point2D> &pointlist,
    bool filterBounds = true )
  {
    pointlist.clear();

    Point2D pt, pt_displaced;
    unsigned importance;
    TileLeafIterator it(root);
    do {
      pt = it->getP1();
      if ( it->isSamplingType() ) // Only "pentagonal" tiles generate sampling points.
      {
        importance = getImportanceAt_bounded( pt );

        // Threshold the function against the F-Code value.
        if ( importance >= calcFCodeValue( it->getFCode() , 2*it->getLevel() ) )
        {
          // Get the displaced point using the lookup table.
          pt_displaced = it->getDisplacedSamplingPoint(importance);

          if ( !filterBounds ||
            (pt_displaced.x>=0 && pt_displaced.x<width &&
            pt_displaced.y>=0 && pt_displaced.y<height) )
          {
            pointlist.push_back(pt_displaced); // collect point.
          }
        }
      }

    } while ( it.next() );

  }

public:

  /// This virtual function must be implemented in order to use the sampling system.
  /// It should return the value of the importance function at the given point.
  virtual unsigned getImportanceAt( Point2D pt ) = 0;

  /// Builds and collects the point set generated be the sampling system,
  /// using the previously defined importance function.
  std::vector<Point2D> getSamplingPoints()
  {
    if (root) delete root;
    std::vector<Point2D> pointlist;

    buildAdaptiveSubdivision();
    collectPoints(pointlist);
    return pointlist;
  }

  /// \example example.cpp
  /// This a simple example of how to use the Quasisampler class.

  /// \example example2.cpp
  /// This another simple example of how to use the Quasisampler class.

  /// \example example3.cpp
  /// This example shows how to use the system on a grayscale image.
};



/*

Static Member initialization

*/

#endif //QUASISAMPLER_PROTOTYPE_H


