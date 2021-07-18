#include <iostream>
#include <fstream>
#include <queue>

using namespace std;

// TODO:
//  - degenerate case when two sites have the same y coord;
//  - degenerate case when circle and site have the same y coord;

struct point {
    double x, y;
};

typedef point site;

struct circle;

struct seg {
    point start, end;
};

struct arc {
    site s;
    circle *c;
    arc *prev, *next;
    seg *l, *r;
};
struct circle{
    double lowest;
    arc *a;
    bool valid;
};

struct compare_events {
    bool operator()(site a, site b) {
        if (a.y < b.y) {
            return true;
        }
        else {
            return false;
        }
    }
    bool operator()(circle *a, circle *b) {
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
    priority_queue<site, vector<site>, compare_events> sites;
    priority_queue<circle*, vector<circle*>, compare_events> circles;


    void handle_site_event() {
        while(!sites.empty()) {
            cout << sites.top().x << " " << sites.top().y << endl;
            sites.pop();
        }
    }

    void handle_circle_event() {
        
    }
    
    public:
        void push(site p) {
            sites.push(p);
        }

        void compute() {

            // While queues are not empty
            while (!sites.empty() || !circles.empty()) {
                if (!sites.empty() && !circles.empty()) {
                    if (sites.top().y > circles.top()->lowest) {
                        handle_site_event();
                    }
                    else {
                        handle_circle_event();
                    }
                }
                else if (!sites.empty()) {
                    handle_site_event();
                }
                else {
                    handle_circle_event();
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
    site p;

    // Initialize site event queue
    while (f >> p.x >> p.y) {
        voronoi.push(p);
    }

    // Compute the Voronoi diagram
    voronoi.compute();

    return 0;
}
