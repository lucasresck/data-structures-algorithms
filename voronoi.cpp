#include <iostream>
#include <fstream>
#include <queue>
#include <math.h> 
#include <random>

using namespace std;

// TODO:
//  - degenerate case when two sites have the same y coord;
//  - degenerate case when circle and site have the same y coord;
//  - handle degeneracy when the point lies exactly bellow intersection between parabolas;
//  - create linked list class;
//  - use Boost's balanced binary search tree;

struct Circle;
struct Vertex;

struct Color {
    int r, g, b;
    
    Color() {};
    Color(int r, int g, int b) : r(r), g(g), b(b) {};
};

/**
 * 2D point, (x, y).
 */

struct Point {
    double x, y;
    Color color {};

    Point() {};
    Point(double x, double y) : x(x), y(y) { };
};

/**
 * Site event for Voronoi algorithm.
 */

typedef Point Site;

/**
 * RGB color.
 */

/**
 * Halfedge of DCEL.
 */

struct Halfedge {
    Vertex *origin;
    Halfedge *prev, *next;
    Halfedge *twin;
    bool cellExtracted;
    Color color;

    Halfedge() : origin(nullptr), prev(nullptr), next(nullptr), twin(nullptr), cellExtracted(false) {};
};

/**
 * Vertex of DCEL.
 */

struct Vertex {
    Point p;
    Halfedge *he;

    Vertex(Point p) : p(p), he(nullptr) {};
};

/**
 * Parabola arc of Voronoi algorithm. Parabola arcs
 * structure is a linked list. In the original
 * implementation, it should be a balanced binary search
 * tree.
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
 * Doubly connected edge list implementation. Store
 * halfedges edges and vertices, which are linked
 * among themselves to form the structure. Faces are
 * not implemented.
 */

class DCEL {
    public:
        vector<Halfedge*> halfedges;
        vector<Vertex*> vertices;

        Halfedge* newHalfedge() {
            Halfedge *he = new Halfedge;
            halfedges.push_back(he);
            return he;
        }

        Vertex* newVertex(Point p) {
            Vertex *v = new Vertex(p);
            vertices.push_back(v);
            return v;
        }

        void linkTwins(Halfedge *he1, Halfedge *he2) {
            he1->twin = he2;
            he2->twin = he1;
        }

        Halfedge* linkVertices(Vertex *v1, Vertex *v2) {
            Halfedge *he1 = newHalfedge();
            he1->origin = v1;
            Halfedge *he2 = newHalfedge();
            he2->origin = v2;
            linkTwins(he1, he2);

            // If vertex 1 has no halfedge.
            if (!(v1->he)) {
                v1->he = he1;
            }

            if (!(v2->he)) {
                v2->he = he2;
            }
            return he1;
        }

        void linkEdges(Halfedge *he1, Halfedge *he2) {
            he1->twin->next = he2;
            he2->prev = he1->twin;

            he2->twin->next = he1;
            he1->prev = he2->twin;
        }

        Point sumPoints(Point p, Point q) {
            Point r;
            r.x = p.x + q.x;
            r.y = p.y + q.y;
            return r;
        }

        Point subPoints(Point p, Point q) {
            Point r;
            r.x = p.x - q.x;
            r.y = p.y - q.y;
            return r;
        }

        double crossPoints(Point v, Point w) {
            return v.x*w.y - v.y*w.x;
        }

        /**
         * Check if two edges intersect. We are considering here that
         * both of them are limited edges, that is, both of them have
         * both extremities. This is usually the case for the first edge
         * (the edge to be cut), because this situation was already tested.
         * For the second edge, it will also be true because this is an edge
         * from the boundign box. Source:
         * https://stackoverflow.com/a/565282/12724988
         * 
         * @param he1 Pointer to halfedge.
         * @param he2 Pointer to halfedge.
         */

        bool checkEdgesIntersection(Halfedge *he1, Halfedge *he2) {
            Point p = he1->origin->p;
            Point r = subPoints(he1->twin->origin->p, p);
            Point q = he2->origin->p;
            Point s = subPoints(he2->twin->origin->p, q);

            // Two colinear or parallel lines
            if (crossPoints(r, s) == 0) return false;

            double t = crossPoints(subPoints(q, p), s)/crossPoints(r, s);
            double u = crossPoints(subPoints(q, p), r)/crossPoints(r, s);

            if ((0 <= t) && (t <= 1) && (0 <= u) && (u <= 1)) return true;

            return false;
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
    DCEL graph;

    // Beginning of linked list
    Arc *arcRoot = nullptr;

    // Bounding box variables
    double bbx1, bbx2, bby1, bby2;
    double margin = 1;

    // Root of bounding box.
    Vertex *boundaryRoot = nullptr;

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
        cout << "Intersection: ";
        if (arc->prev) {
            x1 = calculateIntersection(arc, arc->prev, p.y);
            cout << x1 << " ";
        }
        if (arc->next) {
            x2 = calculateIntersection(arc->next, arc, p.y);
            cout << x2;
        }
        cout << endl;

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
     * Given a parabola arc, check convergence of circle which has
     * the arc as middle arc.
     */

    void checkCircle(Arc *arc) {
        // Check convergence of three breakpoints
        // where arc is the middle arc
        if (arc->prev && arc->next) {
            if (doesConverge(arc->prev->s, arc->s, arc->next->s)) {
                Point center;
                double lowest = calculateLowest(arc->prev->s, arc->s, arc->next->s, &center);
                cout << "Circle event found: " << arc->s.x << " " << lowest << endl;
                Circle *c = new Circle(lowest, arc, center);
                arc->c = c;
                circles.push(c);
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
            Halfedge *he1 = newArc->s0 = newArc->prev->s1 = graph.newHalfedge();
            Halfedge *he2 = newArc->s1 = newArc->next->s0 = graph.newHalfedge();
            graph.linkTwins(he1, he2);

            // Check the "new circles" for convergence.
            checkCircle(newArc->prev);
            checkCircle(newArc->next);
        }
        cout << endl;
    }

    /**
     * Remove an arc of the linked list.
     * 
     * @param a Pointer to arc.
     */

    void removeArc(Arc *a) {
        a->prev->next = a->next;
        a->next->prev = a->prev;
    }

    /**
     * Handle a circle event when this is the next event
     * to be handled.
     */

    void handleCircleEvent() {
        // Get next event and erase it from the queue.
        Circle *c = circles.top();
        circles.pop();

        // There are some circle events in the priority queue
        // which are not valid. We did not remove them because
        // we are dealing with a priority queue, so we just
        // invalidated them.
        if (c->valid) {
            cout << "Handling circle event: ";
            cout << c->a->s.x << " " << c->lowest << endl;

            // Remove the arc from the linked list.
            // However, it still exists.
            removeArc(c->a);

            // Invalidate circle events from successor and predecessor.
            if (c->a->next->c) c->a->next->c->valid = false;
            c->a->next->c = nullptr;
            if (c->a->prev->c) c->a->prev->c->valid = false;
            c->a->prev->c = nullptr;

            // Create vertex at center.
            Vertex *vertex = graph.newVertex(c->center);
            cout << "Vertex: " << c->center.x << " " << c->center.y << endl;
            // Create two halfedge records corresponding to the new
            // breakpoint.
            Halfedge *heDown = graph.newHalfedge();
            Halfedge *heUp = graph.newHalfedge();
            // Set pointers.
            vertex->he = heDown;
            heDown->origin = vertex;
            graph.linkTwins(heDown, heUp);

            // Update breakpoints of left and right arcs.
            c->a->prev->s1 = heUp;
            // We know that heUp must carry the color of the site of the prev of c->a
            heUp->color = c->a->prev->s.color;
            c->a->next->s0 = heDown;
            heDown->color = c->a->next->s.color;

            // Attach the new records to the halfedge records that end
            // at the vertex.

            // Just check which halfedge will be taken from each side,
            // considering the existence of an  origin.
            Halfedge *left = c->a->s0;
            if (c->a->s0->origin) {
                left = c->a->s0->twin;
            }
            left->origin = vertex;
            left->color = c->a->prev->s.color;
            left->twin->color = c->a->s.color;
            Halfedge *right = c->a->s1;
            if (c->a->s1->origin) {
                right = c->a->s1->twin;
            }
            right->origin = vertex;
            right->color = c->a->s.color;
            right->twin->color = c->a->next->s.color;

            // Link among halfedges.
            left->twin->next = right;
            right->prev = left->twin;

            right->twin->next = heDown;
            heDown->prev = right->twin;

            heUp->next = left;
            left->prev = heUp;

            // Check new circles.
            checkCircle(c->a->prev);
            checkCircle(c->a->next);

            delete c->a;
            cout << endl;
        }
        delete c;
    }

    void computeBoundingBox() {
        // Create vertices of bounding box.
        Vertex *a = boundaryRoot = graph.newVertex(Point(bbx1, bby1));
        Vertex *b = graph.newVertex(Point(bbx2, bby1));
        Vertex *d = graph.newVertex(Point(bbx1, bby2));
        Vertex *c = graph.newVertex(Point(bbx2, bby2));

        // Create edges of bounding box.
        Halfedge *ab = graph.linkVertices(a, b);
        Halfedge *bc = graph.linkVertices(b, c);
        Halfedge *cd = graph.linkVertices(c, d);
        Halfedge *da = graph.linkVertices(d, a);


        // Until now, the halfedges of bounding box exist,
        // but they do not have prev and next pointers.
        // We set them now.
        graph.linkEdges(ab, da->twin);
        graph.linkEdges(bc, ab->twin);
        graph.linkEdges(cd, bc->twin);
        graph.linkEdges(da, cd->twin);
    }

    void updateBoundingBox(Point p) {
        double x = p.x;
        double y = p.y;

        if (sites.size() == 1) {
            bbx1 = x - margin;
            bbx2 = x + margin;
            bby1 = y + margin;
            bby2 = y - margin;
        }
        else {
            if (x < bbx1 + margin) {
                bbx1 = x - margin;
            }
            else if (x > bbx2 - margin) {
                bbx2 = x + margin;
            }
            if (y > bby1 - margin) {
                bby1 = y + margin;
            }
            else if (y < bby2 + margin) {
                bby2 = y - margin;
            }
        }
    }

    /** Check if the edge that contains a halfedge
     * has both extremities.
     * 
     * @param he Pointer to halfedge.
     */

    bool isEdgeLimited(Halfedge *he) {
        if (!(he->origin)) return false;
        if (!(he->twin->origin)) return false;
        return true;
    }

    bool isPointBounded(Point p) {
        if (p.x < bbx1 || p.x > bbx2) return false;
        if (p.y > bby1 || p.y < bby2) return false;
        return true;
    }

    /** Check if the (limited) edge that contains a halfedge
     * is inside the bounding box.
     * 
     * @param he Pointer to halfedge.
     * @return the number of bounded extremities. 1 means there is
     * an intersection, for example.
     */

    int isEdgeBounded(Halfedge *he) {
        Point origin1 = he->origin->p;
        Point origin2 = he->twin->origin->p;
        return isPointBounded(origin1) + isPointBounded(origin2);
    }

    /**
     * Intersect an edge that cross the bounding box.
     * 
     * @param he1 Pointer to halfedge.
     * @param he2 Pointer to halfedge of the boundary.
     */

    void intersectEdges(Halfedge *he, Halfedge *heBound) {
    }

    /**
     * Cut an edge that intersects the bounding box.
     * 
     * @param he Pointer to halfedge.
     */

    void cutTheEdge(Halfedge *he) {
        // Iterate over the bounding box.
        Halfedge *heBound = boundaryRoot->he;
        do {
            if (graph.checkEdgesIntersection(he, heBound)) {
                intersectEdges(he, heBound);
                return;
            }
            else {
                heBound = heBound->next;
            }
        }
        while (heBound != boundaryRoot->he);
    }

    void limitDiagramToBoundary() {
        vector<Halfedge*>::iterator he_it;
        for (he_it = graph.halfedges.begin(); he_it != graph.halfedges.end(); he_it++) {
            if (isEdgeLimited(*he_it)) {
                // Number of bounded edge extremities.
                int n = isEdgeBounded(*he_it);
                
                // One point inside, one point outside: intersection.
                if (n == 1) {
                    // It works!
                    cutTheEdge(*he_it);
                }
                // Entire edge outside bounding box. Erase it.
                // Not sure if this can happen.
                else if (n == 0) {
                    //TO DO
                }
            }
            else {
            }
        }        
    }
    
    public:
        /**
         * Push a site event to (site) event queue.
         */
        
        void push(Site p) {
            sites.push(p);
            updateBoundingBox(p);
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

            // computeBoundingBox();
            // limitDiagramToBoundary();

            cout << "Halfedges are: " << endl;
            vector<Halfedge*>::iterator it;
            for (it = graph.halfedges.begin(); it < graph.halfedges.end(); it++) {
                cout << "   ";
                if ((*it)->origin) {
                    cout << (*it)->origin->p.x << " " << (*it)->origin->p.y;
                }
                else {
                    cout << "None" << " " << "None";
                }
                cout << "   ";
                if ((*it)->twin->origin) {
                    cout << (*it)->twin->origin->p.x << " " << (*it)->twin->origin->p.y;
                }
                else {
                    cout << "None" << " " << "None";
                }
                cout << endl;
            }
            cout << endl;

            cout << "Vertices are: " << endl;
            vector<Vertex*>::iterator it2;
            for (it2 = graph.vertices.begin(); it2 != graph.vertices.end(); ++it2) {
                cout << "   " << (*it2)->p.x << " " << (*it2)->p.y << endl;
            }
            cout << endl;

            // Toy examples to test DCEL linking
            // Vertex *vertex = graph.vertices[graph.vertices.size()-1];
            // Point p = vertex->he->next->next->next->next->origin->p;
            // Point p = vertex->he->prev->prev->prev->prev->origin->p;
            // Point p = vertex->he->twin->prev->prev->prev->origin->p;
            // Point p = vertex->he->twin->prev->prev->twin->next->twin->prev->twin->next->next->next->next->next->twin->prev->twin->next->next->next->origin->p;
            // cout << p.x << " " << p.y << endl;
        }

        /**
         * From a specific halfedge, iterate over it until finish
         * the face, and extract the cell that compose this face.
         * 
         * @param begin Pointer to halfedge.
         */

        vector<Point> extractCellFromEdge(Halfedge *begin) {
            // Empty array in case case the face is not well calculated (yet).
            vector<Point> empty {};

            if (begin->cellExtracted) return empty;
            begin->cellExtracted = true;

            // The vector of points.
            if (!(begin->origin)) return empty;
            vector<Point> cell {begin->origin->p};

            // We iterate over the next edges, to extract the points.
            for (Halfedge *he = begin->next; he != begin; he = he->next) {
                if (!he) return empty;
                if (!(he->origin)) return empty;
                Point p = he->origin->p;
                he->cellExtracted = true;
                cell.push_back(p);
            }

            return cell;
        }

        /**
         * Save the extracted cells to file.
         * 
         * @param cells Vector of cells (which are vectors of points).
         */

        void saveCells(vector<vector<Point>> cells, vector<Color> colors) {
            ofstream ofs ("cells.txt", ofstream::out);

            vector<vector<Point>>::iterator cells_it;
            // Iterate over the cells and save them to file.
            int i = 0;
            for (cells_it = cells.begin(); cells_it != cells.end(); cells_it++, i++) {
                vector<Point> cell = *cells_it;
                if (!cell.size()) continue;

                Color color = colors[i];
                ofs << color.r << " " << color.g << " " << color.b << " ";

                vector<Point>::iterator cell_it;
                for (cell_it = cell.begin(); cell_it != cell.end(); cell_it++) {
                    Point p = *cell_it;
                    ofs << p.x << " " << p.y << " ";
                }
                ofs << endl;
            }
            ofs.close();
        }

        /**
         * Extract the cells that form the Voronoi diagram.
         */

        void extractCells() {
            vector<Halfedge*>::iterator he_it;
            vector<vector<Point>> cells;
            vector<Color> colors;
            for (he_it = graph.halfedges.begin(); he_it != graph.halfedges.end(); he_it++) {
                vector<Point> cell = extractCellFromEdge(*he_it);
                cells.push_back(cell);
                colors.push_back((*he_it)->color);
            }
            saveCells(cells, colors);
        }

};

int main() {
    // Random perturbation to avoid degenerate cases.
    mt19937 gen(42);
    uniform_real_distribution<> dis(-0.1, 0.1);

    // Create voronoi object
    Voronoi voronoi;

    // Read the points from the file
    ifstream f("points.txt");

    Site p;
    int r, g, b;

    // Ignore first line, because it has the dimensions
    // of the image.  
    string dummyLine;
    getline(f, dummyLine);

    // Initialize site event queue
    while (f >> r >> g >> b >> p.x >> p.y) {
        // Add little perturbation to avoid degenerate
        // cases when computing Voronoi diagram.
        p.x = p.x + dis(gen);
        p.y = p.y + dis(gen);

        // Save the color of the sites.
        p.color = Color(r, g, b);

        voronoi.push(p);
    }

    // Compute the Voronoi diagram
    voronoi.compute();

    // Extract the cells that form the Voronoi diagram.
    voronoi.extractCells();

    return 0;
}
