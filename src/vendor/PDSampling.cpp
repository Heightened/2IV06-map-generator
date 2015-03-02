#include "PDSampling.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <map>

#include <algorithm>

typedef std::vector<int> IntVector;

#ifndef LARGE_INTEGER
	#include "windows.h"
#endif

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   Modified to be a C++ class by Daniel Dunbar.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* initializes mt[N] with a seed */
RNG::RNG(unsigned long s)
{
	seed(s);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
RNG::RNG(unsigned long init_key[], int key_length)
{
    int i, j, k;
	seed(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

void RNG::seed(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long RNG::getInt32()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long RNG::getInt31()
{
    return (long)(getInt32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double RNG::getDoubleLR()
{
    return getInt32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double RNG::getDoubleL()
{
    return getInt32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double RNG::getDouble()
{
    return (((double)getInt32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

float RNG::getFloatLR()
{
    return getInt32()*(1.0f/4294967295.0f);
    /* divided by 2^32-1 */
}
float RNG::getFloatL()
{
    return getInt32()*(1.0f/4294967296.0f);
    /* divided by 2^32 */
}
float RNG::getFloat()
{
    return (getInt32() + 0.5f)*(1.0f/4294967296.0f);
    /* divided by 2^32 */
}

PDSampler::PDSampler(float _radius, bool _isTiled, bool usesGrid) :
	m_rng((unsigned long) (timeInSeconds()*1000)),
	radius(_radius),
	isTiled(_isTiled)
{
	if (usesGrid) {
			// grid size is chosen so that 4*radius search only
			// requires searching adjacent cells, this also
			// determines max points per cell
		m_gridSize = (int) ceil(2./(4.*_radius));
		if (m_gridSize<2) m_gridSize = 2;

		m_gridCellSize = 2.0f/m_gridSize;
		m_grid = new int[m_gridSize*m_gridSize][kMaxPointsPerCell];

		for (int y=0; y<m_gridSize; y++) {
			for (int x=0; x<m_gridSize; x++) {
				for (int k=0; k<kMaxPointsPerCell; k++) {
					m_grid[y*m_gridSize + x][k] = -1;
				}
			}
		}
	} else {
		m_gridSize = 0;
		m_gridCellSize = 0;
		m_grid = 0;
	}
}

bool PDSampler::pointInDomain(Vec2 &a)
{
	return -1<=a.x && -1<=a.y && 1>=a.x && 1>=a.y;
}

Vec2 PDSampler::randomPoint()
{
	return Vec2(2*m_rng.getFloatL()-1, 2*m_rng.getFloatL()-1);
}

Vec2 PDSampler::getTiled(Vec2 v)
{
	float x = v.x, y = v.y;

	if (isTiled) {
		if (x<-1) x += 2;
		else if (x>1) x -= 2;

		if (y<-1) y += 2;
		else if (y>1) y -= 2;
	}

	return Vec2(x,y);
}

void PDSampler::getGridXY(Vec2 &v, int *gx_out, int *gy_out)
{
	int gx = *gx_out = (int) floor(.5*(v.x + 1)*m_gridSize);
	int gy = *gy_out = (int) floor(.5*(v.y + 1)*m_gridSize);
	if (gx<0 || gx>=m_gridSize || gy<0 || gy>=m_gridSize) {
		printf("Internal error, point outside grid was generated, ignoring.\n");
	}
}

void PDSampler::addPoint(Vec2 pt)
{
	int i, gx, gy, *cell;

	points.push_back(pt);

	if (m_grid) {
		getGridXY(pt, &gx, &gy);
		cell = m_grid[gy*m_gridSize + gx];
		for (i=0; i<kMaxPointsPerCell; i++) {
			if (cell[i]==-1) {
				cell[i] = (int) points.size()-1;
				break;
			}
		}
		if (i==kMaxPointsPerCell) {
			printf("Internal error, overflowed max points per grid cell. Exiting.\n");
			exit(1);
		}
	}
}

int PDSampler::findNeighbors(Vec2 &pt, float distance)
{
	if (!m_grid) {
		printf("Internal error, sampler cannot search without grid.\n");
		exit(1);
	}

	float distanceSqrd = distance*distance;
	int i, j, k, gx, gy, N = (int) ceil(distance/m_gridCellSize);
	if (N>(m_gridSize>>1)) N = m_gridSize>>1;

	m_neighbors.clear();
	getGridXY(pt, &gx, &gy);
	for (j=-N; j<=N; j++) {
		for (i=-N; i<=N; i++) {
			int cx = (gx+i+m_gridSize)%m_gridSize;
			int cy = (gy+j+m_gridSize)%m_gridSize;
			int *cell = m_grid[cy*m_gridSize + cx];

			for (k=0; k<kMaxPointsPerCell; k++) {
				if (cell[k]==-1) {
					break;
				} else {
					if (getDistanceSquared(pt, points[cell[k]])<distanceSqrd)
						m_neighbors.push_back(cell[k]);
				}
			}
		}
	}

	return (int) m_neighbors.size();
}

float PDSampler::findClosestNeighbor(Vec2 &pt, float distance)
{
	if (!m_grid) {
		printf("Internal error, sampler cannot search without grid.\n");
		exit(1);
	}

	float closestSqrd = distance*distance;
	int i, j, k, gx, gy, N = (int) ceil(distance/m_gridCellSize);
	if (N>(m_gridSize>>1)) N = m_gridSize>>1;

	getGridXY(pt, &gx, &gy);
	for (j=-N; j<=N; j++) {
		for (i=-N; i<=N; i++) {
			int cx = (gx+i+m_gridSize)%m_gridSize;
			int cy = (gy+j+m_gridSize)%m_gridSize;
			int *cell = m_grid[cy*m_gridSize + cx];

			for (k=0; k<kMaxPointsPerCell; k++) {
				if (cell[k]==-1) {
					break;
				} else {
					float d = getDistanceSquared(pt, points[cell[k]]);

					if (d<closestSqrd)
						closestSqrd = d;
				}
			}
		}
	}

	return sqrt(closestSqrd);
}

void PDSampler::findNeighborRanges(int index, RangeList &rl)
{
	if (!m_grid) {
		printf("Internal error, sampler cannot search without grid.\n");
		exit(1);
	}

	Vec2 &candidate = points[index];
	float rangeSqrd = 4*4*radius*radius;
	int i, j, k, gx, gy, N = (int) ceil(4*radius/m_gridCellSize);
	if (N>(m_gridSize>>1)) N = m_gridSize>>1;

	getGridXY(candidate, &gx, &gy);

	int xSide = (candidate.x - (-1 + gx*m_gridCellSize))>m_gridCellSize*.5;
	int ySide = (candidate.y - (-1 + gy*m_gridCellSize))>m_gridCellSize*.5;
	int iy = 1;
	for (j=-N; j<=N; j++) {
		int ix = 1;

		if (j==0) iy = ySide;
		else if (j==1) iy = 0;

		for (i=-N; i<=N; i++) {
			if (i==0) ix = xSide;
			else if (i==1) ix = 0;

				// offset to closest cell point
			float dx = candidate.x - (-1 + (gx+i+ix)*m_gridCellSize);
			float dy = candidate.y - (-1 + (gy+j+iy)*m_gridCellSize);

			if (dx*dx+dy*dy<rangeSqrd) {
				int cx = (gx+i+m_gridSize)%m_gridSize;
				int cy = (gy+j+m_gridSize)%m_gridSize;
				int *cell = m_grid[cy*m_gridSize + cx];

				for (k=0; k<kMaxPointsPerCell; k++) {
					if (cell[k]==-1) {
						break;
					} else if (cell[k]!=index) {
						Vec2 &pt = points[cell[k]];
						Vec2 v = getTiled(pt-candidate);
						float distSqrd = v.x*v.x + v.y*v.y;

						if (distSqrd<rangeSqrd) {
							float dist = sqrt(distSqrd);
							float angle = atan2(v.y,v.x);
							float theta = acos(.25f*dist/radius);

							rl.subtract(angle-theta, angle+theta);
						}
					}
				}
			}
		}
	}
}

void PDSampler::maximize()
{
	RangeList rl(0,0);
	int i, N = (int) points.size();

	for (i=0; i<N; i++) {
		Vec2 &candidate = points[i];

		rl.reset(0, (float) M_PI*2);
		findNeighborRanges(i, rl);
		while (rl.numRanges) {
			RangeEntry &re = rl.ranges[m_rng.getInt31()%rl.numRanges];
			float angle = re.min + (re.max-re.min)*m_rng.getFloatL();
			Vec2 pt = getTiled(Vec2(candidate.x + cos(angle)*2*radius,
									candidate.y + sin(angle)*2*radius));

			addPoint(pt);
			rl.subtract(angle - (float) M_PI/3, angle + (float) M_PI/3);
		}
	}
}

void PDSampler::relax()
{
	FILE *tmp = fopen("relaxTmpIn.txt","w");
	int dim, numVerts, numFaces;
	Vec2 *verts = 0;
	int numPoints = (int) points.size();

		// will overwrite later
	fprintf(tmp, "2                  \n");
	for (int i=0; i<(int) points.size(); i++) {
		Vec2 &pt = points[i];
		fprintf(tmp, "%f %f\n", pt.x, pt.y);
	}
	for (int y=-1; y<=1; y++) {
		for (int x=-1; x<=1; x++) {
			if (x || y) {
				for (int i=0; i<(int) points.size(); i++) {
					Vec2 &pt = points[i];
					if (fabs(pt.x+x*2)-1<radius*4 || fabs(pt.y+y*2)-1<radius*4) {
						fprintf(tmp, "%f %f\n", pt.x+x*2, pt.y+y*2);
						numPoints++;
					}
				}
			}
		}
	}
	fseek(tmp, 0, 0);
	fprintf(tmp, "2 %d", numPoints);
	fclose(tmp);

	tmp = fopen("relaxTmpOut.txt", "w");
	fclose(tmp);
	system("qvoronoi p FN < relaxTmpIn.txt > relaxTmpOut.txt");

	tmp = fopen("relaxTmpOut.txt", "r");
	fscanf(tmp, "%d\n%d\n", &dim, &numVerts);

	if (dim!=2) {
		printf("Error calling out to qvoronoi, skipping relaxation.\n");
		goto exit;
	}

	verts = new Vec2[numVerts];
	for (int i=0; i<numVerts; i++) {
		fscanf(tmp, "%f %f\n", &verts[i].x, &verts[i].y);
	}

	fscanf(tmp, "%d\n", &numFaces);

	for (int i=0; i<(int) points.size(); i++) {
		Vec2 center(0,0);
		int N, skip=0;

		fscanf(tmp, "%d", &N);
		for (int j=0; j<N; j++) {
			int index;

			fscanf(tmp, "%d", &index);
			if (index<0) {
				skip = 1;
			} else {
				center += verts[index];
			}
		}

		if (!skip) {
			center *= (1.0f/N);
			points[i] = getTiled(center);
		}
	}

exit:
	if (verts) delete verts;
}

///

DartThrowing::DartThrowing(float radius, bool isTiled, int minMaxThrows, int maxThrowsMult) :
	PDSampler(radius, isTiled),
	m_minMaxThrows(minMaxThrows),
	m_maxThrowsMult(maxThrowsMult)
{
	;
}

void DartThrowing::complete()
{
	while (1) {
		int i, N = (int) points.size()*m_maxThrowsMult;
		if (N<m_minMaxThrows) N = m_minMaxThrows;

		for (i=0; i<N; i++) {
			Vec2 pt = randomPoint();

			findNeighbors(pt, 2*radius);

			if (!m_neighbors.size()) {
				addPoint(pt);
				break;
			}
		}

		if (i==N)
			break;
	}
}


BestCandidate::BestCandidate(float radius, bool isTiled, int multiplier) :
	PDSampler(radius, isTiled),
	m_multiplier(multiplier),
	m_N((int) (.7/(radius*radius)))
{
	;
}

void BestCandidate::complete()
{
	for (int i=0; i<m_N; i++) {
		Vec2 best(0,0);
		float bestDistance = 0;
		int count = 1 + (int) points.size()*m_multiplier;

		for (int j=0; j<count; j++) {
			Vec2 pt = randomPoint();
			float closest = 2;

			closest = findClosestNeighbor(pt, 4*radius);
			if (j==0 || closest>bestDistance) {
				bestDistance = closest;
				best = pt;
			}
		}

		addPoint(best);
	}
}

///

void BoundarySampler::complete()
{
	RangeList rl(0,0);
	IntVector candidates;

	addPoint(randomPoint());
	candidates.push_back((int) points.size()-1);

	while (candidates.size()) {
		int c = m_rng.getInt32()%candidates.size();
		int index = candidates[c];
		Vec2 candidate = points[index];
		candidates[c] = candidates[candidates.size()-1];
		candidates.pop_back();

		rl.reset(0, (float) M_PI*2);
		findNeighborRanges(index, rl);
		while (rl.numRanges) {
			RangeEntry &re = rl.ranges[m_rng.getInt32()%rl.numRanges];
			float angle = re.min + (re.max-re.min)*m_rng.getFloatL();
			Vec2 pt = getTiled(Vec2(candidate.x + cos(angle)*2*radius,
									candidate.y + sin(angle)*2*radius));

			addPoint(pt);
			candidates.push_back((int) points.size()-1);

			rl.subtract(angle - (float) M_PI/3, angle + (float) M_PI/3);
		}
	}
}

///

typedef std::map<int, ScallopedRegion*> RegionMap;

void PureSampler::complete()
{
	Vec2 pt = randomPoint();
	ScallopedRegion *rgn = new ScallopedRegion(pt, radius*2, radius*4);
	RegionMap regions;
	WeightedDiscretePDF<int> regionsPDF;

	addPoint(pt);
	regions[(int) points.size()-1] = rgn;
	regionsPDF.insert((int) points.size()-1, rgn->area);

	while (regions.size()) {
		int idx = regionsPDF.choose(m_rng.getFloatL());

		pt = getTiled(((*regions.find(idx)).second)->sample(m_rng));
		rgn = new ScallopedRegion(pt, radius*2, radius*4);

		findNeighbors(pt, radius*8);
		for (IntVector::const_iterator it=m_neighbors.begin(); it!=m_neighbors.end(); it++) {
			int nIdx = *it;
			Vec2 &n = points[nIdx];

			rgn->subtractDisk(pt+getTiled(n-pt), radius*4);

			RegionMap::iterator entry = regions.find(nIdx);
			if (entry!=regions.end()) {
				ScallopedRegion *nRgn = (*entry).second;
				nRgn->subtractDisk(n+getTiled(pt-n), radius*2);

				if (nRgn->isEmpty()) {
					regions.erase(entry);
					regionsPDF.remove(nIdx);
					delete nRgn;
				} else {
					regionsPDF.update(nIdx, nRgn->area);
				}
			}
		}

		addPoint(pt);

		if (!rgn->isEmpty()) {
			regions[(int) points.size()-1] = rgn;
			regionsPDF.insert((int) points.size()-1, rgn->area);
		} else {
			delete rgn;
		}
	}
}

///

void LinearPureSampler::complete()
{
	IntVector candidates;

	addPoint(randomPoint());
	candidates.push_back((int) points.size()-1);

	while (candidates.size()) {
		int c = m_rng.getInt32()%candidates.size();
		int index = candidates[c];
		Vec2 candidate = points[index];
		candidates[c] = candidates[candidates.size()-1];
		candidates.pop_back();

		ScallopedRegion sr(candidate, radius*2, radius*4);
		findNeighbors(candidate, radius*8);

		for (IntVector::const_iterator it=m_neighbors.begin(); it!=m_neighbors.end(); it++) {
			int nIdx = *it;
			Vec2 &n = points[nIdx];
			Vec2 nClose = candidate + getTiled(n-candidate);

			if (nIdx<index) {
				sr.subtractDisk(nClose, radius*4);
			} else {
				sr.subtractDisk(nClose, radius*2);
			}
		}

		while (!sr.isEmpty()) {
			Vec2 p = sr.sample(m_rng);
			Vec2 pt = getTiled(p);

			addPoint(pt);
			candidates.push_back((int) points.size()-1);

			sr.subtractDisk(p, radius*2);
		}
	}
}

///

class PenroseQuasisampler : public Quasisampler {
	unsigned int val;

public:
	PenroseQuasisampler(unsigned int _val) : Quasisampler(100,100), val(_val) {}

	unsigned int getImportanceAt(Point2D pt) { return val; }
};

void PenroseSampler::complete()
{
	PenroseQuasisampler s((unsigned int) (9.1/(radius*radius)));
	std::vector<Point2D> pts = s.getSamplingPoints();

	for (std::vector<Point2D>::iterator it=pts.begin(); it!=pts.end(); it++ ) {
		Vec2 pt((float) it->x/50.f - 1.0f, (float) it->y/50.f - 1.0f);

		if (pointInDomain(pt)) {
			addPoint(pt);
		}
	}
}

///

void UniformSampler::complete()
{
	int N = (int) (.75/(radius*radius));

	for (int i=0; i<N; i++) {
		addPoint(randomPoint());
	}
}

static const float kTwoPi = (float) (M_PI*2);

static float integralOfDistToCircle(float x, float d, float r, float k)
{
	if (r<FLT_EPSILON)
		return 0.0;

	float sin_x = sin(x);
	float d_sin_x = d*sin_x;
	float y = sin_x*d/r;
	if (y<-1) y = -1;
	else if (y>1) y = 1;

	float theta = asin(y);

	return (r*(r*(x +
				  k*theta) +
			   k*cos(theta)*d_sin_x) +
		    d*cos(x)*d_sin_x)*.5f;
}

ScallopedSector::ScallopedSector(Vec2 &_Pt, float _a1, float _a2, Vec2 &P1, float r1, float sign1, Vec2 &P2, float r2, float sign2)
{
	Vec2 v1 = Vec2(P1.x - _Pt.x, P1.y - _Pt.y);
	Vec2 v2 = Vec2(P2.x - _Pt.x, P2.y - _Pt.y);

	P = _Pt;
	a1 = _a1;
	a2 = _a2;

	arcs[0].P = P1;
	arcs[0].r = r1;
	arcs[0].sign = sign1;
	arcs[0].d = sqrt(v1.x*v1.x + v1.y*v1.y);
	arcs[0].rSqrd = arcs[0].r*arcs[0].r;
	arcs[0].dSqrd = arcs[0].d*arcs[0].d;
	arcs[0].theta = atan2(v1.y,v1.x);
	arcs[0].integralAtStart = integralOfDistToCircle(a1 - arcs[0].theta, arcs[0].d, arcs[0].r, arcs[0].sign);

	arcs[1].P = P2;
	arcs[1].r = r2;
	arcs[1].sign = sign2;
	arcs[1].d = sqrt(v2.x*v2.x + v2.y*v2.y);
	arcs[1].rSqrd = arcs[1].r*arcs[1].r;
	arcs[1].dSqrd = arcs[1].d*arcs[1].d;
	arcs[1].theta = atan2(v2.y,v2.x);
	arcs[1].integralAtStart = integralOfDistToCircle(a1 - arcs[1].theta, arcs[1].d, arcs[1].r, arcs[1].sign);

	area = calcAreaToAngle(a2);
}

float ScallopedSector::calcAreaToAngle(float angle)
{
	float underInner = integralOfDistToCircle(angle - arcs[0].theta, arcs[0].d, arcs[0].r, arcs[0].sign) - arcs[0].integralAtStart;
	float underOuter = integralOfDistToCircle(angle - arcs[1].theta, arcs[1].d, arcs[1].r, arcs[1].sign) - arcs[1].integralAtStart;

	return underOuter-underInner;
}

float ScallopedSector::calcAngleForArea(float area, RNG &rng)
{
	float lo = a1, hi = a2, cur = lo + (hi-lo)*rng.getFloat();

	for (int i=0; i<10; i++) {
		if (calcAreaToAngle(cur)<area) {
			lo = cur;
			cur = (cur + hi)*.5f;
		} else {
			hi = cur;
			cur = (lo + cur)*.5f;
		}
	}

	return cur;
}

float ScallopedSector::distToCurve(float angle, int index)
{
	float alpha = angle - arcs[index].theta;
	float sin_alpha = sin(alpha);
	float t0 = arcs[index].rSqrd - arcs[index].dSqrd*sin_alpha*sin_alpha;
	if (t0<0) {
		return arcs[index].d*cos(alpha);
	} else {
		return arcs[index].d*cos(alpha) + arcs[index].sign*sqrt(t0);
	}
}

Vec2 ScallopedSector::sample(RNG &rng)
{
	float angle = calcAngleForArea(area*rng.getFloatL(), rng);
	float d1 = distToCurve(angle, 0);
	float d2 = distToCurve(angle, 1);
	float d = sqrt(d1*d1 + (d2*d2 - d1*d1)*rng.getFloat());

	return Vec2(P.x + cos(angle)*d, P.y + sin(angle)*d);
}

///

float ScallopedSector::canonizeAngle(float angle)
{
	float delta = fmod(angle - a1, kTwoPi);
	if (delta<0) delta += kTwoPi;
	return a1 + delta;
}

void ScallopedSector::distToCircle(float angle, Vec2 &C, float r, float *d1_out, float *d2_out)
{
	Vec2 v(C.x - P.x, C.y - P.y);
	float dSqrd = v.x*v.x + v.y*v.y;
	float theta = atan2(v.y, v.x);
	float alpha = angle - theta;
	float sin_alpha = sin(alpha);
	float xSqrd = r*r - dSqrd*sin_alpha*sin_alpha;

	if (xSqrd<0) {
		*d1_out = *d2_out = -10000000;
	} else {
		float a = sqrt(dSqrd)*cos(alpha);
		float x = sqrt(xSqrd);
		*d1_out = a-x;
		*d2_out = a+x;
	}
}

void ScallopedSector::subtractDisk(Vec2 &C, float r, std::vector<ScallopedSector> *regions)
{
	std::vector<float> angles;

	Vec2 v(C.x - P.x, C.y-P.y);
	float d = sqrt(v.x*v.x + v.y*v.y);

	if (r<d) {
		float theta = atan2(v.y, v.x);
		float x = sqrt(d*d-r*r);
		float angle, alpha = asin(r/d);

		angle = canonizeAngle(theta+alpha);
		if (a1<angle && angle<a2) {
			if (distToCurve(angle,0)<x && x<distToCurve(angle,1))
				angles.push_back(angle);
		}

		angle = canonizeAngle(theta-alpha);
		if (a1<angle && angle<a2) {
			if (distToCurve(angle,0)<x && x<distToCurve(angle,1))
				angles.push_back(angle);
		}
	}

	for (int arcIndex=0; arcIndex<2; arcIndex++) {
		Vec2 &C2 = arcs[arcIndex].P;
		float R = arcs[arcIndex].r;
		Vec2 v(C.x - C2.x, C.y - C2.y);
		float d = sqrt(v.x*v.x + v.y*v.y);

		if (d>FLT_EPSILON) {
			float invD = 1.0f/d;
			float x = (d*d - r*r + R*R)*invD*.5f;
			float k = R*R - x*x;

			if (k>0) {
				float y = sqrt(k);
				float vx = v.x*invD;
				float vy = v.y*invD;
				float vx_x = vx*x, vy_x = vy*x;
				float vx_y = vx*y, vy_y = vy*y;
				float angle;

				angle = canonizeAngle(atan2(C2.y + vy_x + vx_y - P.y,
											C2.x + vx_x - vy_y - P.x));
				if (a1<angle && angle<a2) angles.push_back(angle);

				angle = canonizeAngle(atan2(C2.y + vy_x - vx_y - P.y,
											C2.x + vx_x + vy_y - P.x));
				if (a1<angle && angle<a2) angles.push_back(angle);
			}
		}
	}

	sort(angles.begin(), angles.end());
	angles.insert(angles.begin(), a1);
	angles.push_back(a2);

	for (unsigned int i=1; i<angles.size(); i++) {
		float a1 = angles[i-1], a2 = angles[i];
		float midA = (a1+a2)*.5f;
		float inner = distToCurve(midA,0);
		float outer = distToCurve(midA,1);
		float d1, d2;

		distToCircle(midA, C, r, &d1, &d2); // d1<=d2

		if (d2<inner || d1>outer) {
			regions->push_back(ScallopedSector(P, a1, a2, arcs[0].P, arcs[0].r, arcs[0].sign, arcs[1].P, arcs[1].r, arcs[1].sign));
		} else {
			if (inner<d1) {
				regions->push_back(ScallopedSector(P, a1, a2, arcs[0].P, arcs[0].r, arcs[0].sign, C, r, -1));
			}
			if (d2<outer) {
				regions->push_back(ScallopedSector(P, a1, a2, C, r, 1, arcs[1].P, arcs[1].r, arcs[1].sign));
			}
		}
	}
}

///

ScallopedRegion::ScallopedRegion(Vec2 &P, float r1, float r2, float _minArea) :
	minArea(_minArea)
{
	regions = new std::vector<ScallopedSector>;
	regions->push_back(ScallopedSector(P, 0, kTwoPi, P, r1, 1, P, r2, 1));
	area = (*regions)[0].area;
}

ScallopedRegion::~ScallopedRegion()
{
	delete regions;
}

void ScallopedRegion::subtractDisk(Vec2 C, float r)
{
	std::vector<ScallopedSector> *newRegions = new std::vector<ScallopedSector>;

	area = 0;
	for (unsigned int i=0; i<regions->size(); i++) {
		ScallopedSector &ss = (*regions)[i];
		std::vector<ScallopedSector> *tmp = new std::vector<ScallopedSector>;

		ss.subtractDisk(C, r, tmp);

		for (unsigned int j=0; j<tmp->size(); j++) {
			ScallopedSector &nss = (*tmp)[j];

			if (nss.area>minArea) {
				area += nss.area;

				if (newRegions->size()) {
					ScallopedSector &last = (*newRegions)[newRegions->size()-1];
					if (last.a2==nss.a1 && (last.arcs[0].P==nss.arcs[0].P && last.arcs[0].r==nss.arcs[0].r && last.arcs[0].sign==nss.arcs[0].sign) &&
						(last.arcs[1].P==nss.arcs[1].P && last.arcs[1].r==nss.arcs[1].r && last.arcs[1].sign==nss.arcs[1].sign)) {
						last.a2 = nss.a2;
						last.area = last.calcAreaToAngle(last.a2);
						continue;
					}
				}

				newRegions->push_back(nss);
			}
		}

		delete tmp;
	}

	delete regions;
	regions = newRegions;
}

Vec2 ScallopedRegion::sample(RNG &rng)
{
	if (!regions->size()) {
		printf("Fatal error, sampled from empty region.");
		exit(1);
		return Vec2(0,0);
	} else {
		float a = area*rng.getFloatL();
		ScallopedSector &ss = (*regions)[0];

		for (unsigned int i=0; i<regions->size(); i++) {
			ss = (*regions)[i];
			if (a<ss.area)
				break;
			a -= ss.area;
		}

		return ss.sample(rng);
	}
}


static const float kSmallestRange = .000001f;

RangeList::RangeList(float min, float max)
{
	numRanges = 0;
	rangesSize = 8;
	ranges = new RangeEntry[rangesSize];
	reset(min, max);
}

RangeList::~RangeList()
{
	delete[] ranges;
}

void RangeList::reset(float min, float max)
{
	numRanges = 1;
	ranges[0].min = min;
	ranges[0].max = max;
}

void RangeList::deleteRange(int pos)
{
	if (pos<numRanges-1) {
		memmove(&ranges[pos], &ranges[pos+1], sizeof(*ranges)*(numRanges-(pos+1)));
	}
	numRanges--;
}

void RangeList::insertRange(int pos, float min, float max)
{
	if (numRanges==rangesSize) {
		RangeEntry *tmp = new RangeEntry[rangesSize];
		memcpy(tmp, ranges, numRanges*sizeof(*tmp));
		delete[] ranges;
		ranges = tmp;
	}

	if (pos<numRanges) {
		memmove(&ranges[pos+1], &ranges[pos], sizeof(*ranges)*(numRanges-pos));
	}

	ranges[pos].min = min;
	ranges[pos].max = max;
	numRanges++;
}

void RangeList::subtract(float a, float b)
{
	static const float twoPi = (float) (M_PI*2);

	if (a>twoPi) {
		subtract(a-twoPi, b-twoPi);
	} else if (b<0) {
		subtract(a+twoPi, b+twoPi);
	} else if (a<0) {
		subtract(0, b);
		subtract(a+twoPi,twoPi);
	} else if (b>twoPi) {
		subtract(a, twoPi);
		subtract(0, b-twoPi);
	} else if (numRanges==0) {
		;
	} else {
		int pos;

		if (a<ranges[0].min) {
			pos = -1;
		} else {
			int lo=0, mid=0, hi=numRanges;

			while (lo<hi-1) {
				mid = (lo+hi)>>1;
				if (ranges[mid].min<a) {
					lo = mid;
				} else {
					hi = mid;
				}
			}

			pos = lo;
		}

		if (pos==-1) {
			pos = 0;
		} else if (a<ranges[pos].max) {
			float c = ranges[pos].min;
			float d = ranges[pos].max;
			if (a-c<kSmallestRange) {
				if (b<d) {
					ranges[pos].min = b;
				} else {
					deleteRange(pos);
				}
			} else {
				ranges[pos].max = a;
				if (b<d) {
					insertRange(pos+1, b, d);
				}
				pos++;
			}
		} else {
			if (pos<numRanges-1 && b>ranges[pos+1].min) {
				pos++;
			} else {
				return;
			}
		}

		while (pos<numRanges && b>ranges[pos].min) {
			if (ranges[pos].max-b<kSmallestRange) {
				deleteRange(pos);
			} else {
				ranges[pos].min = b;
			}
		}
	}
}

void RangeList::print()
{
	printf("[");
	for (int i=0; i<numRanges; i++) {
		printf("(%f,%f)%s", ranges[i].min, ranges[i].max, (i==numRanges-1)?"":", ");
	}
	printf("]\n");
}


#ifdef _WIN32
double timeInSeconds()
{
	static double perfFreq;
	static int hasPerfTimer = -1;

	if (hasPerfTimer==-1) {
		LARGE_INTEGER perfFreqCountsPerSec;

		if (QueryPerformanceFrequency(&perfFreqCountsPerSec)) {
			hasPerfTimer = 1;
			perfFreq = (double) perfFreqCountsPerSec.QuadPart;
		} else {
			hasPerfTimer = 0;
		}
	}

	LARGE_INTEGER count;
	if (hasPerfTimer && QueryPerformanceCounter(&count)) {
		return (double) count.QuadPart/perfFreq;
	} else {
		return GetTickCount()/1000.0;
	}
}
#else
#include <time.h>
#include <sys/time.h>

double timeInSeconds()
{
	struct timeval tv;
	struct timezone tz;

	gettimeofday(&tv, &tz);
	return tv.tv_sec + tv.tv_usec/1000000.0;
}
#endif


const unsigned Quasisampler::fiboTable[32]=
{ 1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,
1597,2584,4181,6765,10946,17711,28657,46368,75025,
121393,196418,317811,514229,832040,1346269,2178309 };

const Point2D Quasisampler::vvect[]={
  Point2D(0,1), Point2D(-0.309017,0.951057), Point2D(-0.587785,0.809017),
  Point2D(-0.809017,0.587785), Point2D(-0.951057,0.309017), Point2D(-1,0),
  Point2D(-0.951057,-0.309017), Point2D(-0.809017,-0.587785),
  Point2D(-0.587785,-0.809017), Point2D(-0.309017,-0.951057), Point2D(0,-1),
  Point2D(0.309017,-0.951057), Point2D(0.587785,-0.809017), Point2D(0.809017,-0.587785),
  Point2D(0.951057,-0.309017), Point2D(1,0), Point2D(0.951057,0.309017),
  Point2D(0.809017,0.587785), Point2D(0.587785,0.809017), Point2D(0.309017,0.951057)
};

const double Quasisampler::lut[LUT_SIZE][21][2] =
{{{0.0130357, 0.0419608}, {-0.0241936, 0.0152706}, {-0.00384601, -0.311212}, {-0.000581893, -0.129134},
  {-0.0363269, 0.0127624}, {0.0999483, 0.408639}, {-0.0526517, 0.4385}, {-0.128703, 0.392}, {0.0132026, 1.0818},
  {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{0.00793289, 0.0148063}, {0.0206067, -0.0809589}, {0.0110103, -0.430433}, {0.0000473169, -0.293185},
  {-0.0593578, 0.019457}, {0.34192, 0.291714}, {-0.286696, 0.386017}, {-0.345313, 0.311961}, {0.00606029, 1.00877},
  {0.04757, 0.05065}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{0.00454493, -0.00805726}, {0.0545058, -0.140953}, {0.00960599, -0.493483}, {0.000527191, -0.354496},
  {-0.0742085, -0.0477178}, {0.436518, 0.218493}, {-0.422435, 0.275524}, {-0.425198, 0.257027},
  {0.0127468, 0.979585}, {0.128363, 0.139522}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}, {0, 0}}, {{-0.0014899, -0.0438403}, {0.122261, -0.229582}, {-0.00497263, -0.580537},
  {-0.00489546, -0.424237}, {-0.107601, -0.133695}, {0.526304, 0.125709}, {-0.558461, 0.0679206},
  {-0.511708, 0.153397}, {0.0271526, 0.950065}, {0.298021, 0.327582}, {-0.00464701, -0.00362132}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0182024, -0.0837012}, {0.226792, -0.318088}, {-0.0416745, -0.663614}, {-0.0253331, -0.455424},
  {-0.159087, -0.20807}, {0.552691, 0.0525824}, {-0.617244, -0.197362}, {-0.561762, 0.00314535},
  {0.0522991, 0.928754}, {0.376689, 0.429912}, {-0.0180693, -0.00792235}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{-0.0308901, -0.108719}, {0.362157, -0.377329},
  {-0.0918077, -0.742776}, {-0.0571567, -0.453854}, {-0.242014, -0.230347}, {0.542952, -0.00542364},
  {-0.614735, -0.35591}, {-0.565238, -0.204834}, {0.084241, 0.900632}, {0.403207, 0.481046},
  {-0.0459391, -0.00743248}, {0.0143212, 0.0776031}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}}, {{-0.0429758, -0.112222}, {0.470514, -0.41007}, {-0.139291, -0.797567}, {-0.0930261, -0.382258},
  {-0.30831, -0.210972}, {0.504387, -0.05265}, {-0.578917, -0.4354}, {-0.545885, -0.40618}, {0.122368, 0.852639},
  {0.377534, 0.476884}, {-0.0712593, 0.0238995}, {0.0349156, 0.248696}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{-0.0297026, -0.0818903}, {0.514634, -0.426843}, {-0.161039, -0.817284},
  {-0.099245, -0.221824}, {-0.359506, -0.135015}, {0.433957, -0.0878639}, {-0.541453, -0.46714},
  {-0.526484, -0.556459}, {0.1735, 0.771396}, {0.353023, 0.455358}, {-0.07854, 0.0885735}, {0.0714601, 0.591673},
  {-0.0147015, 0.0839976}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0204607, -0.0433266}, {0.515056, -0.428386}, {-0.153717, -0.803384}, {-0.0874438, 0.032819},
  {-0.370233, 0.00469937}, {0.331072, -0.0951004}, {-0.507368, -0.487422}, {-0.533403, -0.648977},
  {0.243233, 0.652577}, {0.33663, 0.406983}, {-0.0624495, 0.167064}, {0.0527702, 0.808443}, {-0.0444704, 0.258347},
  {0.030331, -0.00128903}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0184965, 0.00557424}, {0.495666, -0.40889}, {-0.136052, -0.781115}, {-0.0493628, 0.265293},
  {-0.337945, 0.202038}, {0.193353, -0.0835904}, {-0.479971, -0.497456}, {-0.574003, -0.71938},
  {0.32445, 0.514949}, {0.331709, 0.341565}, {-0.034108, 0.244375}, {0.0149632, 0.910353}, {-0.104428, 0.60938},
  {0.0948414, -0.00216379}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0436899, 0.0294207}, {0.469933, -0.372015}, {-0.153852, -0.756531}, {0.00920944, 0.393625},
  {-0.270292, 0.392355}, {0.0540646, -0.0473047}, {-0.466651, -0.492248}, {-0.647575, -0.793479},
  {0.394352, 0.385016}, {0.330852, 0.272582}, {-0.0125759, 0.30811}, {-0.0407447, 0.902855}, {-0.136947, 0.8021},
  {0.227048, -0.0014045}, {0.0261797, 0.0109521}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0602358, 0.0215278}, {0.43301, -0.338538}, {-0.233311, -0.71494}, {0.0916642, 0.433266},
  {-0.173199, 0.474801}, {-0.0384285, 0.024931}, {-0.475596, -0.469989}, {-0.739327, -0.866143},
  {0.440049, 0.277063}, {0.326099, 0.207864}, {-0.00488013, 0.365323}, {-0.0890991, 0.872087},
  {-0.159106, 0.889116}, {0.311406, 0.0126425}, {0.081674, 0.0403966}, {0.01391, 0.00573611}, {0, 0}, {0, 0},
  {0, 0}, {0, 0}, {0, 0}}, {{-0.0723894, -0.00927744}, {0.354855, -0.326512}, {-0.329593, -0.647058},
  {0.169384, 0.42962}, {-0.0250381, 0.472328}, {-0.108748, 0.122704}, {-0.507741, -0.424372},
  {-0.805866, -0.896362}, {0.48306, 0.211626}, {0.314407, 0.142681}, {-0.00348365, 0.415081},
  {-0.125494, 0.836485}, {-0.183247, 0.847226}, {0.366439, 0.0391043}, {0.18978, 0.100287}, {0.0401008, 0.018797},
  {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{-0.0748666, -0.0517059}, {0.237999, -0.333105},
  {-0.391007, -0.558425}, {0.223599, 0.428175}, {0.159284, 0.420084}, {-0.17834, 0.234411}, {-0.553952, -0.353981},
  {-0.821481, -0.848098}, {0.527132, 0.175271}, {0.312397, 0.0908259}, {0.00190795, 0.441568},
  {-0.149358, 0.790424}, {-0.226469, 0.765995}, {0.383259, 0.0740479}, {0.243694, 0.15335}, {0.0901877, 0.0475938},
  {-0.00963625, 0.00819101}, {0, 0}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0862318, -0.0937052}, {0.132383, -0.310846}, {-0.420153, -0.463782}, {0.261956, 0.440763},
  {0.290379, 0.392449}, {-0.264095, 0.349189}, {-0.576491, -0.274722}, {-0.797096, -0.724963},
  {0.565701, 0.153393}, {0.315376, 0.0546255}, {0.0149326, 0.430477}, {-0.167772, 0.702404}, {-0.283244, 0.645617},
  {0.383304, 0.0988087}, {0.248786, 0.17877}, {0.103708, 0.0729573}, {-0.0286781, 0.0298329},
  {-0.00878083, 0.0189161}, {0, 0}, {0, 0}, {0, 0}},
 {{-0.0911025, -0.116785}, {0.058151, -0.268943}, {-0.424486, -0.374671}, {0.288764, 0.470621},
  {0.362681, 0.386055}, {-0.327219, 0.436709}, {-0.585384, -0.202215}, {-0.772145, -0.5936}, {0.580061, 0.135496},
  {0.313963, 0.0305349}, {0.0109925, 0.360967}, {-0.181933, 0.552414}, {-0.300836, 0.508161},
  {0.364265, 0.0976394}, {0.210088, 0.176749}, {0.096516, 0.0958074}, {-0.0658733, 0.0731591},
  {-0.0280071, 0.057776}, {0.0158411, 0.00325704}, {0, 0}, {0, 0}},
 {{-0.0974734, -0.0918732}, {0.0139633, -0.212455}, {-0.406371, -0.282796}, {0.296357, 0.483457},
  {0.381376, 0.39536}, {-0.333854, 0.503081}, {-0.58254, -0.14516}, {-0.763625, -0.49765}, {0.567887, 0.121286},
  {0.30413, 0.0127316}, {-0.00152308, 0.270083}, {-0.191895, 0.352083}, {-0.283727, 0.35145},
  {0.326415, 0.0742237}, {0.163984, 0.15982}, {0.0726181, 0.108651}, {-0.0800514, 0.114725},
  {-0.0673361, 0.138093}, {0.0402953, 0.00961117}, {-0.0193168, 0.0236477}, {0, 0}},
 {{-0.0790912, -0.0163216}, {-0.00448123, -0.162101}, {-0.352873, -0.196134}, {0.271462, 0.449512},
  {0.35836, 0.383875}, {-0.286884, 0.565229}, {-0.550438, -0.0846486}, {-0.75899, -0.42121}, {0.528606, 0.119818},
  {0.280538, 0.00168322}, {-0.0349212, 0.150096}, {-0.171099, 0.193366}, {-0.250974, 0.211407},
  {0.280682, 0.0548899}, {0.126017, 0.143427}, {0.0562988, 0.110436}, {-0.0785227, 0.145239},
  {-0.0937526, 0.190149}, {0.0791086, 0.0227095}, {-0.0545744, 0.0707386}, {0, 0}},
 {{-0.0518157, 0.0510771}, {-0.00760212, -0.128097}, {-0.253754, -0.111841}, {0.205436, 0.354864},
  {0.295866, 0.325402}, {-0.192075, 0.64807}, {-0.4774, -0.00676484}, {-0.722069, -0.332801}, {0.470923, 0.131373},
  {0.244358, -0.00366888}, {-0.0555535, 0.0625726}, {-0.128642, 0.0933316}, {-0.239777, 0.136585},
  {0.234046, 0.0562388}, {0.105223, 0.134278}, {0.0497268, 0.106459}, {-0.0606163, 0.175207},
  {-0.106271, 0.232174}, {0.0538097, 0.0296093}, {-0.122383, 0.16238}, {-0.0113815, 0.0340113}},
 {{-0.0304857, 0.0883196}, {0.00193379, -0.129688}, {-0.148195, -0.0572436}, {0.128477, 0.258454},
  {0.18546, 0.230594}, {-0.120249, 0.694404}, {-0.326488, 0.130702}, {-0.599671, -0.166452}, {0.371228, 0.215584},
  {0.18765, -0.00862734}, {-0.0530754, 0.00501476}, {-0.0781737, 0.0495139}, {-0.215913, 0.0922068},
  {0.202485, 0.0708782}, {0.103985, 0.125369}, {0.0553649, 0.1009}, {-0.0397036, 0.199708}, {-0.0966645, 0.253069},
  {-0.0153489, 0.0350904}, {-0.134291, 0.193388}, {-0.0315258, 0.0780417}},
 {{-0.00909437, 0.0971829}, {0.00766774, -0.145809}, {-0.0755563, -0.0337505}, {0.0700629, 0.188928},
  {0.109764, 0.175155}, {-0.084045, 0.707208}, {-0.200288, 0.246694}, {-0.431284, 0.0136518}, {0.274276, 0.314326},
  {0.138397, -0.0136486}, {-0.033298, -0.019655}, {-0.0429267, 0.0341841}, {-0.195447, 0.0692005},
  {0.188428, 0.0886883}, {0.112392, 0.115937}, {0.0568682, 0.0920568}, {-0.0238131, 0.214855},
  {-0.0754228, 0.259851}, {-0.0881413, 0.0371697}, {-0.127762, 0.194639}, {-0.0700573, 0.173426}}};

