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

    Halfedge() : origin(nullptr), prev(nullptr), next(nullptr), twin(nullptr) {};
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

    Arc(Site s) : s(s), c(nullptr), prev(nullptr), next(nullptr), s0(nullptr), s1(nullptr) {};
};

struct Circle {
    double lowest;
    Arc *a;
    bool valid;

    Circle(double lowest, Arc *a) : valid(true) {};
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
    vector<Halfedge*> halfedges;
    vector<Vertex*> vertices;

    // Beginning of linked list
    Arc *arcRoot = nullptr;

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
        if (arc->prev) { x1 = calculateIntersection(arc, arc->prev, p.y); }
        if (arc->next) { x2 = calculateIntersection(arc->next, arc, p.y); }

        cout << x1 << " " << x2 << endl;
        cout << arc->next << endl;

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
            cout << arcIt << endl;
            if (doesIntersect(arcIt, p)) {
                cout << "Intersect! " << arcIt->s.x << endl;
                return arcIt;
            }
        }
        return nullptr;
    }

    void insertInFrontOf(Arc *a, Arc *b) {
        b->next = a->next;
        b->prev = a;
        if (a->next) { a->next->prev = b; }
        a->next = b;
    }

    Arc* insertArc(Arc *a, Site p) {
        // Copy a to b
        Arc *b = new Arc(a->s);
        b->c = a->c;
        b->prev = a->prev;
        b->next = a->next;
        b->s0 = a->s0;
        b->s1 = a->s1;

        insertInFrontOf(a, b);

        // Create our new arc and insert it
        Arc* newArc = new Arc(p);
        insertInFrontOf(a, newArc);
        return newArc;
    }

    bool doesConverge(Site p, Site q, Site r) {
        // Check left turn of p, q, r
        double D = q.x*r.y + p.x*q.y + p.y*r.x
                - (p.y*q.x + q.y*r.x + r.y*p.x);
        return D < 0;
    }

    double calculateLowest(Site p, Site q, Site r) {
        // Calculate lowest point of a circle
        // containing p, q, and r.
        // Source: https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/

        double  x1 = p.x, y1 = p.y, x2 = q.x,
                y2 = q.y, x3 = r.x, y3 = r.y;

        double x12 = x1 - x2;
        double x13 = x1 - x3;
    
        double y12 = y1 - y2;
        double y13 = y1 - y3;
    
        double y31 = y3 - y1;
        double y21 = y2 - y1;
    
        double x31 = x3 - x1;
        double x21 = x2 - x1;
    
        double sx13 = pow(x1, 2) - pow(x3, 2);    
        double sy13 = pow(y1, 2) - pow(y3, 2);
    
        double sx21 = pow(x2, 2) - pow(x1, 2);
        double sy21 = pow(y2, 2) - pow(y1, 2);
    
        double f = (sx13*x12 + sy13*x12 + sx21*x13 + sy21*x13)
                    /(2*(y31*x12 - y21*x13));
        double g = (sx13*y12 + sy13*y12 + sx21*y13 + sy21*y13)
                    /(2*(x31*y12 - x21*y13));
    
        double c = -pow(x1, 2) - pow(y1, 2) - 2*g*x1 - 2*f*y1;

        double h = -g;
        double k = -f;
        double r2 = h*h + k*k - c;
    
        double radius = sqrt(r2);

        return k - radius;
    }

    void checkCircles(Arc *arc) {
        // Check convergence of three breakpoints
        // where arc is the left breakpoint
        if (arc->next) {
            if (arc->next->next) {
                if (doesConverge(arc->s, arc->next->s, arc->next->next->s)) {
                    double lowest = calculateLowest(arc->s, arc->next->s, arc->next->next->s);
                    Circle *c = new Circle(lowest, arc->next);
                    arc->next->c = c;
                    circles.push(c);
                }
            }
        }

        // Check right
        if (arc->prev) {
            if (arc->prev->prev) {
                if (doesConverge(arc->prev->prev->s, arc->prev->s, arc->s)) {
                    double lowest = calculateLowest(arc->prev->prev->s, arc->prev->s, arc->s);
                    Circle *c = new Circle(lowest, arc->prev);
                    arc->prev->c = c;
                    circles.push(c);
                }
            }
        }
    }

    void handleSiteEvent() {
        cout << "Handling site event" << endl;
        // Get next event and erase it
        Site p = sites.top();
        sites.pop();

        // If the linked list is empty
        if (arcRoot == nullptr) {
            cout << "Empty linked list" << endl;
            arcRoot = new Arc(p);
        }
        else {
            Arc *a = findIntersection(p);

            cout << a->c << endl;

            // If the arc has a circle event
            // it is a false alarm
            if (a->c) a->c->valid = false;
            a->c = nullptr;

            // Insert the new arc of site p into the beach line,
            // splitting a into two elements in the linked list.
            Arc* newArc = insertArc(a, p);

            // Create new half edge records for the edge separating
            // arc a and arc with site p, traced by the two new breakpoints.
            Halfedge *he1 = newArc->s0 = newArc->prev->s1 = new Halfedge;
            Halfedge *he2 = newArc->s1 = newArc->next->s0 = new Halfedge;
            he1->twin = he2;
            he2->twin = he1;
            halfedges.push_back(he1); halfedges.push_back(he2);

            // Check the "new circles" for convergence.
            checkCircles(newArc);
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
