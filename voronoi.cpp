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

    Arc(Site s) : s(s), c(NULL), prev(NULL), next(NULL), s0(NULL), s1(NULL) {};
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

    void insertArc(Arc *a, Site p) {
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

            insertArc(a, p);
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
