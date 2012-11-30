//
//  main.m
//  Test
//
//  Created by Robert-Jan on 11/29/12.
//  Copyright (c) 2012 CrossProduct. All rights reserved.
//

#import <Foundation/Foundation.h>

#include <vector>
#include <algorithm>
#include <limits>

#define REAL double

namespace delaunay2D
{
    const int DELAUNAY2D_UNDEFINED = -1;
    
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
        
        bool operator <(Point &pointB)
        {
            return (this->x < pointB.x) || (this->x == pointB.x && this->y < pointB.y);
        }

        REAL Distance(Point &pointB)
        {
            REAL deltaX = pointB.x - this->x;
            REAL deltaY = pointB.y - this->y;
            
            return sqrt((deltaX * deltaX) + (deltaY * deltaY));
        }

        REAL DistanceSquared(Point &pointB)
        {
            REAL deltaX = pointB.x - this->x;
            REAL deltaY = pointB.y - this->y;
            
            return ((deltaX * deltaX) + (deltaY * deltaY));
        }
        
        static REAL CrossProduct(Point &p1, Point &p2, Point &p3)
        {
            REAL u1 =  p2.x - p1.x;
            REAL v1 =  p2.y - p1.y;
            REAL u2 =  p3.x - p1.x;
            REAL v2 =  p3.y - p1.y;
            
            return u1 * v2 - v1 * u2;
        }
    };

    struct Edge
    {
        int s, t;
        int l, r;
    
        Edge()
        {
            s = t = l = r = 0;
        }
        
        Edge(int s, int t)
        {
            this->s = s;
            this->t = t;
        }

        Edge(int s, int t, int l, int r)
        {
            this->s = s;
            this->t = t;
            this->l = l;
            this->r = r;
        }
    };
    
    struct Triangle
    {
        int a, b, c;
        
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
        
        void CircumCircle(Point &p1, Point &p2, Point &p3)
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
    };

    struct Delaunay
    {
        std::vector<Point> points;
        std::vector<Edge> edges;
        std::vector<Triangle> triangles;
        
        Circle bC;

        void Triangulate()
        {
            int s = 0, t = 0;
            
            // Initialise.
            int nFaces = 0;
                        
            // Find closest neighbours and add edge to triangulation.
            this->FindClosestNeighbours(this->points, s, t);
            
            // Create seed edge and add it to the triangulation.
            this->AddEdge(s, t, DELAUNAY2D_UNDEFINED, DELAUNAY2D_UNDEFINED);

            int currentEdge = 0;
            while (currentEdge < this->edges.size())
            {
                Triangle triangle;
                
                if (this->edges[currentEdge].l == DELAUNAY2D_UNDEFINED)
                {
                    if (this->CompleteFacet(currentEdge, nFaces, triangle))
                    {
                        this->triangles.push_back(triangle);
                    }
                    
//                    animate(triangulationState);
                }
                
                if (this->edges[currentEdge].r == DELAUNAY2D_UNDEFINED)
                {
                    if (this->CompleteFacet(currentEdge, nFaces, triangle))
                    {
                        this->triangles.push_back(triangle);
                    }
                    
//                    animate(triangulationState);
                }
                
                currentEdge++;
            }
        }

        void FindClosestNeighbours(std::vector<Point> points, int &u, int &v)
        {
            int s = 0, t = 0;
            
            REAL min = std::numeric_limits<REAL>::max();
            
            for (int i = 0; i < points.size() - 1; i++)
            {
                for (int j = i + 1; j < points.size(); j++)
                {
                    REAL d = points[i].DistanceSquared(points[j]);
                    
                    if (d < min)
                    {
                        s = i;
                        t = j;
                        
                        min = d;
                    }
                }
                
            }

            u = s;
            v = t;
        }
        
        int FindEdge(int s, int t)
        {
            bool edgeExists = false;
            
            int i;
            
            for (i = 0; i < this->edges.size(); i++)
            {
                if ((this->edges[i].s == s && this->edges[i].t == t) || (this->edges[i].s == t && this->edges[i].t == s))
                {
                    edgeExists = true;
                    
                    break;
                }
            }
            
            if (edgeExists)
            {
                return i;
            }
            else
            {
                return DELAUNAY2D_UNDEFINED;
            }
        }
        
        int AddEdge(int s, int t, int l, int r)
        {
            int e = this->FindEdge(s, t);
            
            if (e == DELAUNAY2D_UNDEFINED)
            {
                if (s < t)
                {
                    this->edges.push_back(Edge(s, t, l, r));
                    
                    return (int)this->edges.size() - 1;
                }
                else
                {
                    this->edges.push_back(Edge(t, s, r, l));
                    
                    return (int)this->edges.size() - 1;
                }
            }
            else
            {
                return DELAUNAY2D_UNDEFINED;
            }
        }
        
        /*
         * Update the left face of an edge.
         */
        void UpdateLeftFace(int eI, int s, int t, int f)
        {
            if (!((this->edges[eI].s == s && this->edges[eI].t == t) ||
                  (this->edges[eI].s == t && this->edges[eI].t == s)))
            {
                exit(-1);
            //    Panic.panic("updateLeftFace: adj. matrix and edge table mismatch");
            }
            
            if (this->edges[eI].s == s && this->edges[eI].l == DELAUNAY2D_UNDEFINED)
            {
                this->edges[eI].l = f;
            }
            else if (this->edges[eI].t == s && this->edges[eI].r == DELAUNAY2D_UNDEFINED)
            {
                this->edges[eI].r = f;
            }
            else
            {
                exit(-1);
            }
                //Panic.panic("updateLeftFace: attempt to overwrite edge info");
        }
        
        /*
         * Complete a facet by looking for the circle free point to the left
         * of the edge "e_i".  Add the facet to the triangulation.
         *
         * This function is a bit long and may be better split.
         */
        bool CompleteFacet(int eI, int nFaces, Triangle &triangle)
        {
            //    Edge e[] = tri.edge;
        //    RealPoint p[] = tri.point;
            
            
            int s = 0, t = 0;
            
            // Cache s and t.
            if (this->edges[eI].l == DELAUNAY2D_UNDEFINED)
            {
                s = this->edges[eI].s;
                t = this->edges[eI].t;
            }
            else if (this->edges[eI].r == DELAUNAY2D_UNDEFINED)
            {
                s = this->edges[eI].t;
                t = this->edges[eI].s;
            }
            else
            {
                return false;
            }
                        
            // Find a point on left of edge.
            int u;
            
            for (u = 0; u < this->points.size(); u++)
            {
                if (u == s || u == t)
                {
                    continue;
                }
                
                if (Point::CrossProduct(this->points[s], this->points[t], this->points[u]) > 0.0)
                {
                    break;
                }
            }
            
            // Find best point on left of edge.
            int bP = u;
            
            if (bP < this->points.size())
            {
                this->bC.CircumCircle(this->points[s], this->points[t], this->points[bP]);
                
//                animate(triangleState);
                
                for (u = bP + 1; u < this->points.size(); u++)
                {
                    if (u == s || u == t)
                    {
                        continue;
                    }
                    
//                    animate(pointState);
                    
                    REAL cP = Point::CrossProduct(this->points[s], this->points[t], this->points[u]);
                    
                    if (cP > 0.0)
                    {
                        if (bC.Inside(this->points[u]))
                        {
//                            animate(insideState);
                            bP = u;
                            
                            bC.CircumCircle(this->points[s], this->points[t], this->points[u]);
                            
//                            animate(triangleState);
                        }
                    }
                }
            }
            
            // Add new triangle or update edge info if s-t is on hull.
            if (bP < this->points.size())
            {
                // Update face information of edge being completed. 
                this->UpdateLeftFace(eI, s, t, nFaces);
                
                nFaces++;
                
                // Add new edge or update face info of old edge. 
                eI = this->FindEdge(bP, s);
                
                if (eI == DELAUNAY2D_UNDEFINED)
                {
                    // New edge.
                    eI = this->AddEdge(bP, s, nFaces, DELAUNAY2D_UNDEFINED);
                }
                else
                {
                    // Old edge.
                    this->UpdateLeftFace(eI, bP, s, nFaces);
                }
                
                // Add new edge or update face info of old edge. 
                eI = this->FindEdge(t, bP);
                
                if (eI == DELAUNAY2D_UNDEFINED)
                {
                    // New edge.
                    eI = this->AddEdge(t, bP, nFaces, DELAUNAY2D_UNDEFINED);
                }
                else
                {
                    // Old edge.
                    this->UpdateLeftFace(eI, t, bP, nFaces);
                }
            }
            else
            {
                this->UpdateLeftFace(eI, s, t, 0);
            }
            
            triangle = Triangle(s, t, bP);
            
            return true;
        }
    };
}

int main(int argc, const char * argv[])
{
    @autoreleasepool
    {
        NSLog(@"Hello, World!");
        
        delaunay2D::Delaunay delaunay;
        
        delaunay.points.push_back(delaunay2D::Point(-10, -10));
        delaunay.points.push_back(delaunay2D::Point( 10, -10));
        delaunay.points.push_back(delaunay2D::Point(-10,  10));
        delaunay.points.push_back(delaunay2D::Point( 10,  10));
        delaunay.points.push_back(delaunay2D::Point(  3,   2));
        
        delaunay.Triangulate();
        
        for (int i = 0; i < delaunay.edges.size(); i++)
        {
            NSLog(@"%i %i %i %i", delaunay.edges[i].s, delaunay.edges[i].t, delaunay.edges[i].r, delaunay.edges[i].l);
        }
        
        for (int i = 0; i < delaunay.triangles.size(); i++)
        {
            NSLog(@"Triangle: %i %i %i", delaunay.triangles[i].a, delaunay.triangles[i].b, delaunay.triangles[i].c);
        }
        
//        delaunay.points.push
    }
    
    return 0;
}

