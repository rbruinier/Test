//
// modification of source of : http://www.codeguru.com/cpp/cpp/algorithms/general/article.php/c8901/Delaunay-Triangles.htm
// algorithm is the same
//

#import <Foundation/Foundation.h>

#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <string>

#define REAL double

namespace delaunay2D
{
    struct Point
    {
        REAL x;
        REAL y;
        
        Point()
        {
            this->x = 0.0;
            this->y = 0.0;
        }

        Point(REAL x, REAL y)
        {
            this->x = x;
            this->y = y;
        }
        
        friend bool operator < (const Point &pointA, const Point &pointB)
        {
            return (pointA.x < pointB.x) || (pointA.x == pointB.x && pointA.y < pointB.y);
        }
        
        friend bool operator == (const Point &pointA, const Point &pointB)
        {
            return (pointA.x == pointB.x) && (pointA.y == pointB.y);
        }

        REAL Distance(const Point &pointB) const
        {
            REAL deltaX = pointB.x - this->x;
            REAL deltaY = pointB.y - this->y;
            
            return sqrt((deltaX * deltaX) + (deltaY * deltaY));
        }

        REAL DistanceSquared(const Point &pointB) const
        {
            REAL deltaX = pointB.x - this->x;
            REAL deltaY = pointB.y - this->y;
            
            return ((deltaX * deltaX) + (deltaY * deltaY));
        }
        
        static REAL CrossProduct(const Point &p1, const Point &p2, const Point &p3)
        {
            REAL u1 =  p2.x - p1.x;
            REAL v1 =  p2.y - p1.y;
            REAL u2 =  p3.x - p1.x;
            REAL v2 =  p3.y - p1.y;
            
            return u1 * v2 - v1 * u2;
        }
        
        void Description() const
        {
            std::cout << "x = " << this->x << ", y = " << this->y;
            
            std::cout << "\n";
        }
    };

    struct Edge
    {
        const Point *p1;
        const Point *p2;
            
        Edge(const Point *p1, const Point *p2)
        {
            this->p1 = *p1 < *p2 ? p1 : p2;
            this->p2 = *p1 < *p2 ? p2 : p1;
        }
        
        friend bool operator < (const Edge &edgeA, const Edge &edgeB)
        {
            if (*edgeA.p1 == *edgeB.p1)
            {
                return *edgeA.p2 < *edgeB.p2;
            }
                
            return *edgeA.p1 < *edgeB.p1;
        }
                
        friend bool operator == (const Edge &edgeA, const Edge &edgeB)
        {
            return (edgeA.p1 == edgeB.p1) && (edgeA.p2 == edgeB.p2);
        }
    };
        
    struct Circle
    {
        Point center;
        REAL radius;
        
        Circle()
        {
            this->center = Point();
            this->radius = 0.0;
        }
        
        Circle(Point center, REAL radius)
        {
            this->center = center;
            this->radius = radius;
        }
        
        void Set(Point center, REAL radius)
        {
            this->center = center;
            this->radius = radius;
        }
        
        bool Inside(Point &point)
        {
            return (this->center.DistanceSquared(point) < (this->radius * this->radius));
        }
        
        void CircumCircle(const Point &p1, const Point &p2, const Point &p3)
        {
            REAL cp = Point::CrossProduct(p1, p2, p3);
            
            if (cp != 0.0)
            {
                REAL p1Sq = p1.x * p1.x + p1.y * p1.y;
                REAL p2Sq = p2.x * p2.x + p2.y * p2.y;
                REAL p3Sq = p3.x * p3.x + p3.y * p3.y;
                
                REAL num1 = p1Sq * (p2.y - p3.y) + p2Sq * (p3.y - p1.y) + p3Sq * (p1.y - p2.y);
                REAL cx = num1 / (2.0 * cp);
                
                REAL num2 = p1Sq * (p3.x - p2.x) + p2Sq * (p1.x - p3.x) + p3Sq * (p2.x - p1.x);
                REAL cy = num2 / (2.0 * cp);
                
                this->center = Point(cx, cy);
            }
            
            // Radius 
            this->radius = center.Distance(p1);
        }
        
        void Description()
        {
            std::cout << "center = ";
            
            this->center.Description();
            
            std::cout << ", radius = " << this->radius << "\n";
        }
    };

    
    struct Triangle
    {
        int a, b, c;
        
        const Point *points[3];
        
        Circle circle;
        
        Triangle()
        {
            this->a = this->b = this->c;
        }
        
        Triangle(int a, int b, int c)
        {
            this->a = a;
            this->b = b;
            this->c = c;
        }
        
        Triangle(const Point *p1, const Point *p2, const Point *p3)
        {
            this->points[0] = p1;
            this->points[1] = p2;
            this->points[2] = p3;
            
            circle.CircumCircle(*p1, *p2, *p3);
        }
        
        friend bool operator < (const Triangle &triA, const Triangle &triB)
        {
            if (triA.circle.center.x == triB.circle.center.x)
            {
                return triA.circle.center.y < triB.circle.center.y;
            }
            
            return triA.circle.center.x < triB.circle.center.x;
        }
        
        bool IsLeftOf(const Point &vIt) const
        {
            return vIt.x > (this->circle.center.x + this->circle.radius);
        }
        
        void Description() const
        {
            std::cout << "Triangle BEGIN:\n";
            
            std::cout << "Point 1 = " << this->points[0]->x << ", " << this->points[0]->y << "\n";
            std::cout << "Point 2 = " << this->points[1]->x << ", " << this->points[1]->y << "\n";
            std::cout << "Point 3 = " << this->points[2]->x << ", " << this->points[2]->y << "\n";
            
            std::cout << "Circle center = " << this->circle.center.x << ", " << this->circle.center.y << "; radius = " << this->circle.radius << "\n";

            std::cout << "Triangle END\n";
        }
        
        bool Encompasses(const Point &vIt) const
        {
            return vIt.DistanceSquared(this->circle.center) <= this->circle.radius * this->circle.radius;
        }
        
        bool Contains(const Point &vIt) const
        {
            return (this->points[0] == &vIt) || (this->points[1] == &vIt) || (this->points[2] == &vIt);
        }
    };
    

    struct Delaunay
    {
    public:
        std::set<Point> points;
        std::vector<Triangle> triangles;
        
        Delaunay()
        {
            this->superTriangle = NULL;
        }
        
        void Triangulate()
        {
            if (points.size() < 3)
            {
                return;
            }

            this->CreateSuperTriangle();
            
            std::set<Triangle> workset;
            
            workset.insert(*this->superTriangle);
                        
            for (std::set<Point>::iterator vIt = this->points.begin(); vIt != this->points.end(); vIt++)
            {
                // remove all completed triangles from workset
                this->FinishCompletedTrianglesFromWorkset(*vIt, workset);
                
                // find all triangles
                std::set<Edge> removedEdges;
                
                // find all triangles that are affected by this vertex
                std::set<Triangle>::iterator worksetIterator = workset.begin();
                
                while (worksetIterator != workset.end())
                {
                    std::set<Triangle>::iterator next = worksetIterator;
                    
                    next++;

                    if (worksetIterator->Encompasses(*vIt))
                    {
                        
                        workset.erase(worksetIterator);
                        
                        this->HandleEdge(Edge(worksetIterator->points[0], worksetIterator->points[1]), removedEdges);
                        this->HandleEdge(Edge(worksetIterator->points[1], worksetIterator->points[2]), removedEdges);
                        this->HandleEdge(Edge(worksetIterator->points[2], worksetIterator->points[0]), removedEdges);
                    }
                    
                    worksetIterator = next;
                }
                
                // create new triangles from the edges and the current vertex.
                for (std::set<Edge>::iterator eIt = removedEdges.begin(); eIt != removedEdges.end(); eIt++)
                {
                    Triangle triangle = Triangle(eIt->p1, eIt->p2, &(*vIt));
                    
                    workset.insert(triangle);
                }
            }
            
            for (std::set<Triangle>::iterator wIt = workset.begin(); wIt != workset.end(); wIt++)
            {
                if (!wIt->Contains(*this->superTriangle->points[0]) && !wIt->Contains(*this->superTriangle->points[1]) && !wIt->Contains(*this->superTriangle->points[2]))
                {
                    this->triangles.push_back(*wIt);
                }
            }
        }
        
    private:
        Triangle *superTriangle;

        void HandleEdge(const Edge &edge, std::set<Edge> &edges)
        {
            // try to find an equal edge.. if found, remove it from the list because it is a double edge, otherwise add it to the list
            std::set<Edge>::iterator found = edges.find(edge);
            
            if (found == edges.end())
            {
                edges.insert(edge);
            }
            else
            {
                edges.erase(found);
            }
        }
        
        // walks through all points to create bounding box, then performs some calculations to make a super triangle encompassing all points
        void CreateSuperTriangle()
        {
            std::set<Point>::const_iterator vIt = this->points.begin();
            
            REAL xMin = vIt->x;
            REAL yMin = vIt->y;
            
            REAL xMax = xMin;
            REAL yMax = yMin;
            
            vIt++;
            
            for (; vIt != this->points.end(); vIt++)
            {
                xMax = vIt->x;
                
                REAL y = vIt->y;
                
                if (y < yMin) yMin = y;
                if (y > yMax) yMax = y;
            }
            
            REAL dx = xMax - xMin;
            REAL dy = yMax - yMin;
            
            REAL ddx = dx * 0.01;
            REAL ddy = dy * 0.01;
            
            xMin -= ddx;
            xMax += ddx;
            
            dx += 2 * ddx;
            
            yMin -= ddy;
            yMax += ddy;
            
            dy += 2 * ddy;
            
            const REAL SQRT3 = 1.732050808;

            Point *p1 = new Point(xMin - dy * SQRT3 / 3.0, yMin);
            Point *p2 = new Point(xMax + dy * SQRT3 / 3.0, yMin);
            Point *p3 = new Point((xMin + xMax) * 0.5, yMax + dx * SQRT3 * 0.5);
            
            this->superTriangle = new Triangle(p1, p2, p3);
        }
        
        void FinishCompletedTrianglesFromWorkset(const Point &point, std::set<Triangle> &workset)
        {
            std::set<Triangle>::iterator wIt = workset.begin();
            
            while (wIt != workset.end())
            {
                std::set<Triangle>::iterator next = wIt;
                
                next++;
                
                if (wIt->IsLeftOf(point))
                {
                    workset.erase(wIt);
                    
                    if (!wIt->Contains(*this->superTriangle->points[0]) && !wIt->Contains(*this->superTriangle->points[1]) && !wIt->Contains(*this->superTriangle->points[2]))
                    {
                        this->triangles.push_back(*wIt);
                    }
                    
                }
                
                wIt = next;
            }
            
        }
    };
}

int main(int argc, const char * argv[])
{
    @autoreleasepool
    {
        NSLog(@"Hello, World!");
        
        delaunay2D::Delaunay delaunay;

        const int nrOfPoints = 1000;
        
        for (int i = 0; i < nrOfPoints; i++)
        {
            delaunay.points.insert(delaunay2D::Point(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX));
        }
        
        NSLog(@"Start!");
        
//        delaunay.points.insert(delaunay2D::Point(-10, -10));
//        delaunay.points.insert(delaunay2D::Point( 10, -10));
//        delaunay.points.insert(delaunay2D::Point(-10,  10));
//////        delaunay.points.push_back(delaunay2D::Point( 10,  10));
//        delaunay.points.insert(delaunay2D::Point(  3,   2));
//
        delaunay.Triangulate();
        
        NSLog(@"Done: %i triangles", (int)delaunay.triangles.size());
    }
    
    return 0;
}

