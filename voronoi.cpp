#include <iostream>
#include <fstream>
#include <queue>

using namespace std;

// TODO:
//  - degenerate case when two sites have the same y coord;
//  - degenerate case when circle and site have the same y coord;

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
    

    void handleSiteEvent() {
        while(!sites.empty()) {
            cout << sites.top().x << " " << sites.top().y << endl;
            sites.pop();
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
