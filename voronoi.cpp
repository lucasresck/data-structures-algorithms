#include <iostream>
#include <fstream>
#include <queue>
#include <math.h> 

using namespace std;

// TODO:
//  - degenerate case when two sites have the same y coord;
//  - degenerate case when circle and site have the same y coord;
//  - handle degeneracy when the point lies exactly bellow intersection between parabolas;

struct Point {
    double x, y;
};

typedef Point Site;

struct Circle;
struct Vertex;

struct Halfedge {
    Vertex *origin;
    Halfedge *prev, *next;
    Halfedge *twin;
};

struct Vertex {
    Point p;
    Halfedge *he;
};

struct Arc {
    Site s;
    Circle *c;
    Arc *prev, *next;
    Halfedge *s0, *s1;
};

struct Circle {
    double lowest;
    Arc *a;
    bool valid;
};

struct CompareEvents {
    bool operator()(Site a, Site b) {
        if (a.y < b.y) {
            return true;
        }
        else {
            return false;
        }
    }
    bool operator()(Circle *a, Circle *b) {
        if (a->lowest < b->lowest) {
            return true;
        }
        else {
            return false;
        }
    }
};

class Voronoi {
    // Event queues
    // Instead of only one queue for both events, we are using
    // two, one for each "type" of event. Because of this,
    // we will only need to compare the first element of each queue.
    priority_queue<Site, vector<Site>, CompareEvents> sites;
    priority_queue<Circle*, vector<Circle*>, CompareEvents> circles;

    // Beginning of linked list
    Arc *arcRoot;

    double calculateIntersection(Arc *arc0, Arc *arc1, double l) {
        // Capture (x, y) for each parabola arc
        Point focus1 = arc1->s;
        double x1 = focus1.x;
        double y1 = focus1.y;
        Point focus0 = arc0->s;
        double x0 = focus0.x;
        double y0 = focus0.y;

        double c0 = y0 - l;
        double c1 = y1 - l;

        double a = c1 - c0;
        double b = -2*(x0*c1 - x1*c0);
        double c = x0*x0*c1 + y0*y0*c1 - x1*x1*c0 - y1*y1*c0 - l*l*(c1 - c0);

        return (-b - sqrt(b*b - 4*a*c))/(2*a);
        // Recall that parabolas can intersect at two points.
        // However, considering only -b - sqrt(), and not
        // -b + sqrt(), is enough for our bussiness.
    }

    bool doesIntersect(Arc *arc, Point p) {
        // Calculate x intersection of arc with prev and next
        double x1, x2;
        if (arc->prev) { x1 = calculateIntersection(arc->prev, arc, p.y); }
        if (arc->next) { x2 = calculateIntersection(arc, arc->next, p.y); }

        // Verify if p.x is between both intersections
        // That is, if p intercept arc
        // If there is not arc->prev, for example,
        // we will have x1 = -infinity
        if ((!(arc->prev) || x1 <= p.x) && (!(arc->next) || p.x <= x2)) { return true; }
        else { return false; }
    }

    Arc* findIntersection(Site p) {
        // Iterate over the beach line
        for (Arc *arcIt = arcRoot; arcIt; arcIt = arcIt->next) {
            if (doesIntersect(arcIt, p)) {
                return arcIt;
            }
        }
    }

    void handleSiteEvent() {
        // Get next event and erase it
        Site p = sites.top();
        sites.pop();

        // If the linked list is empty
        if (arcRoot == nullptr) {
            arcRoot = new Arc;
            arcRoot->s = p;
        }
        else {
            Arc *a = findIntersection(p);
        }
    }

    void handleCircleEvent() {
        
    }
    
    public:
        void push(Site p) {
            sites.push(p);
        }

        void compute() {

            // While queues are not empty
            while (!sites.empty() || !circles.empty()) {
                if (!sites.empty() && !circles.empty()) {
                    if (sites.top().y > circles.top()->lowest) {
                        handleSiteEvent();
                    }
                    else {
                        handleCircleEvent();
                    }
                }
                else if (!sites.empty()) {
                    handleSiteEvent();
                }
                else {
                    handleCircleEvent();
                }
            }
        }

};

int main() {
    // Utilization example

    // Create voronoi object
    Voronoi voronoi;

    // Read the points from the file
    ifstream f("points.txt");
    Site p;

    // Initialize site event queue
    while (f >> p.x >> p.y) {
        voronoi.push(p);
    }

    // Compute the Voronoi diagram
    voronoi.compute();

    return 0;
}
