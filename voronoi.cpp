#include <iostream>
#include <fstream>
#include <queue>
#include <math.h> 

using namespace std;

// TODO:
//  - degenerate case when two sites have the same y coord;
//  - degenerate case when circle and site have the same y coord;
//  - handle degeneracy when the point lies exactly bellow intersection between parabolas;

/**
 * 2D point, (x, y).
 */

struct Point {
    double x, y;
};

/**
 * Site event for Voronoi algorithm.
 */

typedef Point Site;

struct Circle;
struct Vertex;

/**
 * Halfedge of DCEL.
 */

struct Halfedge {
    Vertex *origin;
    Halfedge *prev, *next;
    Halfedge *twin;

    Halfedge() : origin(nullptr), prev(nullptr), next(nullptr), twin(nullptr) {};
};

/**
 * Vertex of DCEL.
 */

struct Vertex {
    Point p;
    Halfedge *he;
};

/**
 * Parabola arc of Voronoi algorithm.
 */

struct Arc {
    Site s;
    Circle *c;
    Arc *prev, *next;
    Halfedge *s0, *s1;

    Arc(Site s) : s(s), c(nullptr), prev(nullptr), next(nullptr), s0(nullptr), s1(nullptr) {};
};

/**
 * Circle event of Voronoi algorithm.
 */

struct Circle {
    double lowest;
    Arc *a;
    bool valid;
    Point center;

    Circle(double lowest, Arc *a, Point center) :
    lowest(lowest), a(a), valid(true), center(center) {};
};

/**
 * Comparison of two Voronoi events by y coordinate.
 */

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

/**
 * Voronoi algorithm implementation.
 * Specifically, Fortune's algorithm.
 */

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

    /**
     * Calculate x coordinate of intersection between
     * two parabola arcs.
     * 
     * @param arc0 Pointer to arc.
     * @param arc1 Pointer to arc.
     * @param l y coordinate of sweep line.
     */

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

        // Recall that parabolas can intersect at two points.
        // However, considering only -b - sqrt(), and not
        // -b + sqrt(), is enough for our bussiness.
        return (-b - sqrt(b*b - 4*a*c))/(2*a);
    }

    /**
     * Check if parabola arc is just above point p. We need to find
     * which arc is above p, and this function check one
     * specific arc.
     * 
     * @param arc Pointer to arc.
     * @param p Point p.
     */

    bool doesIntersect(Arc *arc, Point p) {
        // Calculate x intersection of arc with prev and next
        double x1, x2;
        if (arc->prev) { x1 = calculateIntersection(arc, arc->prev, p.y); }
        if (arc->next) { x2 = calculateIntersection(arc->next, arc, p.y); }

        cout << "Intersection: " << x1 << " " << x2 << endl;

        // Verify if p.x is between both intersections
        // That is, if p intercept arc
        // If there is not arc->prev, for example,
        // we will have x1 = -infinity
        if ((!(arc->prev) || x1 <= p.x) && (!(arc->next) || p.x <= x2)) { return true; }
        else { return false; }
    }

    /**
     * Find the parabola arc that is just above site event p.
     * 
     * @param p Site event.
     */

    Arc* findIntersection(Site p) {
        // Iterate over the beach line
        for (Arc *arcIt = arcRoot; arcIt; arcIt = arcIt->next) {
            if (doesIntersect(arcIt, p)) {
                cout << "Intersect with arc with site: " << arcIt->s.x << " " << arcIt->s.y << endl;
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

    /**
     * Create arc for site p and insert it between
     * arc a and a copy of it.
     * 
     * @param a Pointer to arc a.
     * @param p Site event p.
     */

    Arc* insertArc(Arc *a, Site p) {
        // Copy a to b
        Arc *b = new Arc(a->s);
        // b->c = a->c;
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

    /**
     * Check convergence of three points.
     * In Voronoi algorithm, three points p, q, and r
     * can generate a circle event, where the parabola arc
     * related to q will disappear. However, it is possible
     * to check first if this is possible. If it is possible,
     * not necessarily will happen, but it saves computation.
     */

    bool doesConverge(Site p, Site q, Site r) {
        // Check left turn of p, q, r
        double D = q.x*r.y + p.x*q.y + p.y*r.x
                - (p.y*q.x + q.y*r.x + r.y*p.x);
        return D < 0;
    }

    /**
     * Calculate lowest point of circle that intersects points
     * p, q, and r. This is the point we are interested during
     * the algorithm, because it indicates the "priority" of
     * the circle event among all site and circle events.
     * Source: https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/
     * 
     * @param p First point.
     * @param q Second point.
     * @param r Third point.
     * @param center Pointer to center of the arc. This value does
     * not exist yet, in fact it will be "returned".
     */

    double calculateLowest(Site p, Site q, Site r, Point *center) {
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

        center->x = h;
        center->y = k;

        return k - radius;
    }

    /**
     * Given a parabola arc, check convergence of neighbor circles.
     * Check convergence of circles generated by (arc, arc->next,
     * arc->next->next) and (aarc->prev->prev, arc->prev, arc).
     */

    void checkCircles(Arc *arc) {
        // Check convergence of three breakpoints
        // where arc is the left breakpoint
        if (arc->next) {
            if (arc->next->next) {
                if (doesConverge(arc->s, arc->next->s, arc->next->next->s)) {
                    Point center;
                    double lowest = calculateLowest(arc->s, arc->next->s, arc->next->next->s, &center);
                    cout << "Circle event: " << arc->next->s.x << " " << lowest << endl;
                    Circle *c = new Circle(lowest, arc->next, center);
                    arc->next->c = c;
                    circles.push(c);
                }
            }
        }

        // Check right
        if (arc->prev) {
            if (arc->prev->prev) {
                if (doesConverge(arc->prev->prev->s, arc->prev->s, arc->s)) {
                    Point center;
                    double lowest = calculateLowest(arc->prev->prev->s, arc->prev->s, arc->s, &center);
                    cout << "Circle event found: " << arc->prev->s.x << " " << lowest << endl;
                    Circle *c = new Circle(lowest, arc->prev, center);
                    arc->prev->c = c;
                    circles.push(c);
                }
            }
        }
    }

    /**
     * Handle a site event when this is the next event
     * to be handled.
     */

    void handleSiteEvent() {
        cout << "Handling site event: ";
        // Get next event and erase it
        Site p = sites.top();
        sites.pop();
        cout << p.x << " " << p.y << endl;

        // If the linked list is empty
        if (arcRoot == nullptr) {
            cout << "Empty linked list. Adding first arc." << endl;
            arcRoot = new Arc(p);
        }
        else {
            Arc *a = findIntersection(p);

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
        cout << endl;
    }

    /**
     * Delete arc of linked list.
     * 
     * @param a Pointer to arc.
     * @param aNext Pointer to pointer to next arc. This is necessary
     * to continue tracking original arc neighbors after deleting arc a.
     */

    void deleteArc(Arc *a, Arc **aNext) {
        a->prev->next = a->next;
        a->next->prev = a->prev;
        *aNext = a->next;
        delete a;
    }

    /**
     * Handle a circle event when this is the next event
     * to be handled.
     */

    void handleCircleEvent() {
        // Get next event and erase it
        Circle *c = circles.top();
        circles.pop();

        // There are some circle events in the priority queue
        // which are not valid. We did not remove them because
        // we are dealing with a priority queue, so we just
        // invalidated them.
        if (c->valid) {
            cout << "Handling circle event: ";
            cout << c->a->s.x << " " << c->lowest << endl;

            // Save pointer to next arc, because this one
            // will be deleted.
            Arc *aNext;
            deleteArc(c->a, &aNext);

            // Invalidate circle events from successor and predecessor.
            if (aNext->c) aNext->c->valid = false;
            aNext->c = nullptr;
            if (aNext->prev->c) aNext->prev->c->valid = false;
            aNext->prev->c = nullptr;
        }
    }
    
    public:
        /**
         * Push a site event to (site) event queue.
         */
        
        void push(Site p) {
            sites.push(p);
        }

        /**
         * Compute whole Voronoi algorithm.
         */

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
