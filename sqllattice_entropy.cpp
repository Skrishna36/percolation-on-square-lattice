#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <random>
#include <chrono>
#include <numeric>
#include <cmath>

using namespace std;

const double L = 1.0;
const double EPS_COORD = 1e-9;

struct ROW {
    int id;
    // tuple: (y_min, y_max, x_min, x_max)
    vector<tuple<double, double, double, double>> points;
    vector<int> up_neighbors;    // neighbor ids (0-based)
    vector<int> down_neighbors;
    vector<int> left_neighbors;
    vector<int> right_neighbors;
    int area;
};

// Container to return step-by-step tracking results back to main()
struct PercolationResult {
    int spanning_cluster_size;
    vector<double> entropy;
    vector<double> Smax;
};

// thread-local RNG used throughout (avoids RAND_MAX / reseed issues)
static thread_local std::mt19937_64 rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());
static inline double my_random() {
    uniform_real_distribution<double> d(0.0, 1.0);
    return d(rng);
}

// Check if two cells "touch" (share a boundary)
int neighbor_check(const vector<ROW>& rows, int p, int l) {
    if (p < 0 || p >= (int)rows.size() || l < 0 || l >= (int)rows.size()) return 0;
    const ROW& row_p = rows[p];
    const ROW& row_l = rows[l];

    double a = get<0>(row_p.points[0]);
    double b = get<1>(row_p.points[0]);
    double c = get<2>(row_p.points[0]);
    double d = get<3>(row_p.points[0]);

    double a1 = get<0>(row_l.points[0]);
    double b1 = get<1>(row_l.points[0]);
    double c1 = get<2>(row_l.points[0]);
    double d1 = get<3>(row_l.points[0]);

    // Check strict overlap on orthogonal axis and touching on one axis
    bool x_overlap = (c < d1) && (c1 < d);
    bool y_overlap = (a < b1) && (a1 < b);

    if ((fabs(b - a1) < EPS_COORD && x_overlap) ||
        (fabs(a - b1) < EPS_COORD && x_overlap) ||
        (fabs(c - d1) < EPS_COORD && y_overlap) ||
        (fabs(d - c1) < EPS_COORD && y_overlap)) {
        return 1;
    }
    return 0;
}

int findroot(int x, vector<int>& parent) {
    if (parent[x] < 0) return x;
    return parent[x] = findroot(parent[x], parent);
}

void unite(int a, int b,
           vector<int>& parent, vector<int>& redefined_cluster,
           vector<int>& boundary_mask)
{
    a = findroot(a, parent);
    b = findroot(b, parent);
    if (a == b) return;

    // union by size (negative values store size)
    if (parent[a] < parent[b]) {  // a larger
        parent[a] += parent[b];
        redefined_cluster[a] += redefined_cluster[b];
        redefined_cluster[b] = 0;
        parent[b] = a;
        boundary_mask[a] |= boundary_mask[b];
    } else {
        redefined_cluster[b] += redefined_cluster[a];
        redefined_cluster[a] = 0;
        parent[b] += parent[a];
        parent[a] = b;
        boundary_mask[b] |= boundary_mask[a];
    }
}

PercolationResult percolator(int total_cells, const vector<ROW>& rows)
{
    ofstream tout("percolation_n.txt");
    ofstream chosenFile("chosen_cells.txt");

    vector<int> parent(total_cells, 0);
    vector<double> entropy(total_cells + 1, 0.0);
    vector<double> Smax(total_cells + 1, 0.0);
    vector<int> redefined_cluster(total_cells, 0);
    vector<int> boundary_mask(total_cells, 0);
    vector<bool> occupied(total_cells, false);
    vector<vector<int>> nn(total_cells);

    const int TOUCH_TOP    = 1 << 0;
    const int TOUCH_BOTTOM = 1 << 1;
    const int TOUCH_LEFT   = 1 << 2;
    const int TOUCH_RIGHT  = 1 << 3;

    double total_bonds = 0;

    for (int i = 0; i < total_cells; ++i) {
        const ROW& cell = rows[i];
        vector<int> neighbors;

        for (int n : cell.up_neighbors)    neighbors.push_back(n);
        for (int n : cell.down_neighbors)  neighbors.push_back(n);
        for (int n : cell.left_neighbors)  neighbors.push_back(n);
        for (int n : cell.right_neighbors) neighbors.push_back(n);

        nn[i] = neighbors;
        total_bonds += neighbors.size();
        
        tout << i << "\t";
        for (size_t j = 0; j < nn[i].size(); ++j) {
            tout << nn[i][j];
            if (j != nn[i].size() - 1) tout << ",";
        }
        tout << "\n";
    }
    total_bonds = total_bonds / 2.0;

    vector<int> order(total_cells);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);

    int spanning_cluster_size = 0;

    // OCCUPATION LOOP
    for (int t = 0; t < total_cells; ++t)
    {
        int s = order[t];
        occupied[s] = true;
        parent[s] = -1;

        chosenFile << rows[s].id << endl;
        for (const auto& point : rows[s].points) {
            double y_min, y_max, x_min, x_max;
            tie(y_min, y_max, x_min, x_max) = point;
            chosenFile << x_min << " " << y_min << "\n"
                       << x_max << " " << y_min << "\n"
                       << x_max << " " << y_max << "\n"
                       << x_min << " " << y_max << "\n\n";
        }

        // Boundary detection
        auto& p = rows[s].points[0];
        double y_min = get<0>(p), y_max = get<1>(p);
        double x_min = get<2>(p), x_max = get<3>(p);

        if (fabs(y_min - 0.0) < EPS_COORD) boundary_mask[s] |= TOUCH_BOTTOM;
        if (fabs(y_max - L) < EPS_COORD)   boundary_mask[s] |= TOUCH_TOP;
        if (fabs(x_min - 0.0) < EPS_COORD) boundary_mask[s] |= TOUCH_LEFT;
        if (fabs(x_max - L) < EPS_COORD)   boundary_mask[s] |= TOUCH_RIGHT;

        // UNION WITH NEIGHBORS
        int cluster_add = 0;
        redefined_cluster[s] = nn[s].size();
        
        for (int nb : nn[s]) {
            if (occupied[nb]) {
                cluster_add += 1;
                int r1 = findroot(s, parent);
                int r2 = findroot(nb, parent);
                if (r1 != r2) {
                    unite(s, nb, parent, redefined_cluster, boundary_mask);
                }
            }
        }
        int r1 = findroot(s, parent);
        redefined_cluster[r1] = redefined_cluster[r1] - cluster_add;

        // CHECK SPANNING
        int root = findroot(s, parent);
        int mask = boundary_mask[root];
        bool horizontal = (mask & TOUCH_LEFT) && (mask & TOUCH_RIGHT);
        bool vertical   = (mask & TOUCH_TOP) && (mask & TOUCH_BOTTOM);

        if ((horizontal || vertical) && spanning_cluster_size == 0) {
            spanning_cluster_size = t + 1; 
        }

        // ENTROPY AND SMAX CALCULATION FOR STEP (t+1)
        double entropy_step = 0;
        double s_max = 0;
        double unclustered_bonds = total_bonds;

        for (int k = 0; k < total_cells; k++) {
            if (redefined_cluster[k] > 0) {
                if (redefined_cluster[k] > s_max) {
                    s_max = redefined_cluster[k];
                }
                double S = static_cast<double>(redefined_cluster[k]) / total_bonds;
                unclustered_bonds -= redefined_cluster[k]; 
                entropy_step -= (S * log(S));
            }
        }  

        entropy[t + 1] = entropy_step; 
        if (unclustered_bonds > 0) {
            entropy[t + 1] -= unclustered_bonds * ((1.0 / total_bonds) * log(1.0 / total_bonds));
        }
        Smax[t + 1] = s_max;                                
    }

    tout.close();
    chosenFile.close();

    // Pack inside result struct and return to main()
    return {spanning_cluster_size, entropy, Smax};
}

int cellno(double x, double y, const vector<ROW>& rows) {
    int j = -1;
    for (int i = 0; i < (int)rows.size(); i++) {
        const auto& bounds = rows[i].points[0];
        double y_min = get<0>(bounds), y_max = get<1>(bounds);
        double x_min = get<2>(bounds), x_max = get<3>(bounds);
        
        if ((x_min - EPS_COORD < x) && (x < x_max + EPS_COORD) &&
            (y_min - EPS_COORD < y) && (y < y_max + EPS_COORD)) {
            j = i;
            break;
        }
    }
    return j;
}

vector<ROW> lattice(int cut, int total_cells, vector<int> nums) {
    vector<ROW> rows;
    rows.reserve(total_cells + 1);
    rows.push_back({0, {make_tuple(0.0, L, 0.0, L)}, {}, {}, {}, {}, static_cast<int>(L * L)});

    for (int e = 0; e < cut; ++e) {
        int selected_index = cellno(my_random(), my_random(), rows);
        ROW saved_cell = rows[selected_index];
        int selected_id = saved_cell.id;

        auto old_point = saved_cell.points[0];
        double y_min = get<0>(old_point), y_max = get<1>(old_point);
        double x_min = get<2>(old_point), x_max = get<3>(old_point);

        double cut_x = x_min + (x_max - x_min) * my_random();
        double cut_y = y_min + (y_max - y_min) * my_random();

        for (int up_id : saved_cell.up_neighbors) {
            auto& dn = rows[up_id].down_neighbors;
            dn.erase(remove(dn.begin(), dn.end(), selected_id), dn.end());
        }
        for (int down_id : saved_cell.down_neighbors) {
            auto& upn = rows[down_id].up_neighbors;
            upn.erase(remove(upn.begin(), upn.end(), selected_id), upn.end());
        }
        for (int left_id : saved_cell.left_neighbors) {
            auto& right_nbrs = rows[left_id].right_neighbors;
            right_nbrs.erase(remove(right_nbrs.begin(), right_nbrs.end(), selected_id), right_nbrs.end());
        }
        for (int right_id : saved_cell.right_neighbors) {
            auto& left_nbrs = rows[right_id].left_neighbors;
            left_nbrs.erase(remove(left_nbrs.begin(), left_nbrs.end(), selected_id), left_nbrs.end());
        }

        if (nums[e] == 3) {
            int new_id = rows.size();
            double area_bottom_left = (cut_x - x_min) * (cut_y - y_min);
            double area_bottom_right = (x_max - cut_x) * (cut_y - y_min);
            double area_top_left = (cut_x - x_min) * (y_max - cut_y);
            double area_top_right = (x_max - cut_x) * (y_max - cut_y);

            rows[selected_index].points[0] = make_tuple(y_min, cut_y, x_min, cut_x);
            rows[selected_index].area = (area_bottom_left);
            rows[selected_index].right_neighbors = {new_id};  
            rows[selected_index].up_neighbors =  {new_id + 2};     
            rows[selected_index].left_neighbors.clear();
            rows[selected_index].down_neighbors.clear();

            ROW bottom_right;
            bottom_right.id = {new_id};
            bottom_right.points = {make_tuple(y_min, cut_y, cut_x, x_max)};
            bottom_right.area = (area_bottom_right);
            bottom_right.left_neighbors = {selected_id};       
            bottom_right.up_neighbors =  {new_id + 1};            
            bottom_right.right_neighbors.clear();
            bottom_right.down_neighbors.clear();

            ROW top_right;
            top_right.id = {new_id + 1};
            top_right.points = {make_tuple(cut_y, y_max, cut_x, x_max)};
            top_right.area = (area_top_right);
            top_right.left_neighbors =  {new_id + 2};             
            top_right.down_neighbors =  {new_id};             
            top_right.right_neighbors.clear();
            top_right.up_neighbors.clear();

            ROW top_left;
            top_left.id = {new_id + 2};
            top_left.points = {make_tuple(cut_y, y_max, x_min, cut_x)};
            top_left.area = static_cast<int>(area_top_left);
            top_left.right_neighbors =  {new_id + 1};             
            top_left.down_neighbors = {selected_id};          
            top_left.left_neighbors.clear();
            top_left.up_neighbors.clear();

            rows.push_back(bottom_right);
            rows.push_back(top_right);
            rows.push_back(top_left);

            vector<int> new_cells = {selected_id, {new_id}, {new_id + 1}, {new_id + 2}};
            
            for (int cell_index : new_cells) {
                for (int left_id : saved_cell.left_neighbors) {
                    if (neighbor_check(rows, cell_index, left_id)) {
                        rows[cell_index].left_neighbors.push_back(left_id);
                        rows[left_id].right_neighbors.push_back(cell_index);
                    }
                }
                for (int right_id : saved_cell.right_neighbors) {
                    if (neighbor_check(rows, cell_index, right_id)) {
                        rows[cell_index].right_neighbors.push_back(right_id);
                        rows[right_id].left_neighbors.push_back(cell_index);
                    }
                }
                for (int up_id : saved_cell.up_neighbors) {
                    if (neighbor_check(rows, cell_index, up_id)) {
                        rows[cell_index].up_neighbors.push_back(up_id);
                        rows[up_id].down_neighbors.push_back(cell_index);
                    }
                }
                for (int down_id : saved_cell.down_neighbors) {
                    if (neighbor_check(rows, cell_index, down_id)) {
                        rows[cell_index].down_neighbors.push_back(down_id);
                        rows[down_id].up_neighbors.push_back(cell_index);
                    }
                }
            }
        }
        else {
            bool vertical_cut = (my_random() < 0.5);
            if (vertical_cut) {
                double area_left  = (y_max - y_min) * (cut_x - x_min);
                double area_right = (y_max - y_min) * (x_max - cut_x);
                int new_id = rows.size();

                rows[selected_index].points[0] = make_tuple(y_min, y_max, x_min, cut_x);
                rows[selected_index].area = area_left;
                rows[selected_index].left_neighbors.clear();
                rows[selected_index].right_neighbors = {new_id};
                rows[selected_index].up_neighbors.clear();
                rows[selected_index].down_neighbors.clear();

                ROW right_cell;
                right_cell.id = new_id;
                right_cell.points = {make_tuple(y_min, y_max, cut_x, x_max)};
                right_cell.area = area_right;
                right_cell.left_neighbors = {selected_id};
                rows.push_back(right_cell);

                vector<int> new_cells = {selected_id, new_id};
                for (int cell_index : new_cells) {
                    for (int left_id : saved_cell.left_neighbors) {
                        if (neighbor_check(rows, cell_index, left_id)) {
                            rows[cell_index].left_neighbors.push_back(left_id);
                            rows[left_id].right_neighbors.push_back(cell_index);
                        }
                    }
                    for (int right_id : saved_cell.right_neighbors) {
                        if (neighbor_check(rows, cell_index, right_id)) {
                            rows[cell_index].right_neighbors.push_back(right_id);
                            rows[right_id].left_neighbors.push_back(cell_index);
                        }
                    }
                    for (int up_id : saved_cell.up_neighbors) {
                        if (neighbor_check(rows, cell_index, up_id)) {
                            rows[cell_index].up_neighbors.push_back(up_id);
                            rows[up_id].down_neighbors.push_back(cell_index);
                        }
                    }
                    for (int down_id : saved_cell.down_neighbors) {
                        if (neighbor_check(rows, cell_index, down_id)) {
                            rows[cell_index].down_neighbors.push_back(down_id);
                            rows[down_id].up_neighbors.push_back(cell_index);
                        }
                    }
                }
            } else {
                double area_bottom = (x_max - x_min) * (cut_y - y_min);
                double area_top    = (x_max - x_min) * (y_max - cut_y);
                int new_id = rows.size();

                rows[selected_index].points[0] = make_tuple(y_min, cut_y, x_min, x_max);
                rows[selected_index].area = area_bottom;
                rows[selected_index].down_neighbors.clear();
                rows[selected_index].up_neighbors = {new_id};
                rows[selected_index].left_neighbors.clear();
                rows[selected_index].right_neighbors.clear();

                ROW top_cell;
                top_cell.id = new_id;
                top_cell.points = {make_tuple(cut_y, y_max, x_min, x_max)};
                top_cell.area = area_top;
                top_cell.down_neighbors = {selected_id};
                rows.push_back(top_cell);

                vector<int> new_cells = {selected_id, new_id};
                for (int cell_index : new_cells) {
                    for (int left_id : saved_cell.left_neighbors) {
                        if (neighbor_check(rows, cell_index, left_id)) {
                            rows[cell_index].left_neighbors.push_back(left_id);
                            rows[left_id].right_neighbors.push_back(cell_index);
                        }
                    }
                    for (int right_id : saved_cell.right_neighbors) {
                        if (neighbor_check(rows, cell_index, right_id)) {
                            rows[cell_index].right_neighbors.push_back(right_id);
                            rows[right_id].left_neighbors.push_back(cell_index);
                        }
                    }
                    for (int up_id : saved_cell.up_neighbors) {
                        if (neighbor_check(rows, cell_index, up_id)) {
                            rows[cell_index].up_neighbors.push_back(up_id);
                            rows[up_id].down_neighbors.push_back(cell_index);
                        }
                    }
                    for (int down_id : saved_cell.down_neighbors) {
                        if (neighbor_check(rows, cell_index, down_id)) {
                            rows[cell_index].down_neighbors.push_back(down_id);
                            rows[down_id].up_neighbors.push_back(cell_index);
                        }
                    }
                }
            }
        }
    }
    return rows;
}

int totalcells(int cut, vector<int>& nums) {
    int total = 1;
    for (int e = 0; e < cut; ++e) { 
        if (my_random() < 0.5) {
            nums.push_back(1);
            total += 1;
        } else {
            nums.push_back(3);
            total += 3;
        }
    }
    return total;
}

int main() {
    int cuts = 10000;
    int iterations = 100; // Increased to show the effect of averaging
    int maximum_cell = cuts * 3 + 1; 
    
    // Define a uniform continuous probability spectrum for statistical accumulation
    int num_bins = maximum_cell + 3;
    vector<double> probability_spectrum;
    vector<double> percolation_stats(num_bins, 0.0); 
    vector<double> global_avg_entropy(num_bins, 0.0);
    vector<double> global_avg_smax(num_bins, 0.0);
    vector<int> bin_counts(num_bins, 0); // Tracks how many samples fell into each probability bin
 vector<int> nums;
        int total_cells = totalcells(cuts, nums); 
        vector<ROW> rows = lattice(cuts, total_cells, nums);  

    for (int j = 0; j < num_bins; ++j) {
        probability_spectrum.push_back(static_cast<double>(j) / (num_bins - 1));
    }

    for (int trial = 1; trial <= iterations; ++trial) {
       
        
        // Execute percolation step tracking and capture returned results 
        PercolationResult res = percolator(total_cells, rows);
        int cells_to_span = res.spanning_cluster_size;
        double spanning_prob = static_cast<double>(cells_to_span) / total_cells;

        // 1. Accumulate spanning probability statistics
        if (cells_to_span > 0 && cells_to_span <= total_cells) {
            for (size_t i = 0; i < probability_spectrum.size(); ++i) {
                if (probability_spectrum[i] >= spanning_prob) {
                    percolation_stats[i] += 1.0;
                }
            }
        }

        // 2. Recompute total bonds to correctly scale the order parameter output
        double total_bonds = 0;
        for(int i = 0; i < total_cells; ++i) {
            total_bonds += (rows[i].up_neighbors.size() + rows[i].down_neighbors.size() + 
                            rows[i].left_neighbors.size() + rows[i].right_neighbors.size());
        }
        total_bonds /= 2.0;

        // 3. Map discrete microcanonical values to the continuous global bins
        for (int f = 1; f <= total_cells; f++) {
            double occupation_fraction = static_cast<double>(f) / total_cells;
            double normalized_entropy = res.entropy[f] / res.entropy[1];
            double normalized_smax = res.Smax[f] / total_bonds;

            // Find the closest index in the global probability spectrum
            int bin_idx = round(occupation_fraction * (num_bins - 1));
            if (bin_idx >= 0 && bin_idx < num_bins) {
                global_avg_entropy[bin_idx] += normalized_entropy;
                global_avg_smax[bin_idx] += normalized_smax;
                bin_counts[bin_idx]++;
            }
        }

        if (trial % max(1, iterations / 10) == 0) {
            cout << "Completed " << trial << " / " << iterations << endl;
        }
    }

    // --- ORCHESTRATE GLOBAL TRIAL AVERAGING AND WRITE TO FILES ---

    // Finalize the averages by dividing by the accumulated bin counts
    for (int i = 0; i < num_bins; ++i) {
        if (bin_counts[i] > 0) {
            global_avg_entropy[i] /= bin_counts[i];
            global_avg_smax[i] /= bin_counts[i];
        }
    }

    // Write averaged entropy and order parameter (Smax) to entropy.txt
    ofstream eout("entropy.txt");
    if (eout.is_open()) {
        eout << "OccupationFraction AveragedEntropy AveragedSmax\n";
        for (int i = 0; i < num_bins; ++i) {
            // Only output points that were visited during the trials
            if (bin_counts[i] > 0) {
                eout << probability_spectrum[i] << "    " 
                     << global_avg_entropy[i] << "    " 
                     << global_avg_smax[i] << "\n";
            }
        }
        eout.close();
    }

    // Write spanning probability frequency curve to p8000.txt
    ofstream fout("p8000.txt");
    fout << "FractionOccupied Frequency" << endl;
    for (size_t i = 1; i < probability_spectrum.size(); ++i) {
        fout << probability_spectrum[i] << " " 
             << (percolation_stats[i] / iterations) << "\n";
    }
    fout.close();

    cout << "\nPercolation Statistics Summary:" << endl;
    cout << "Total trials: " << iterations << endl;
    cout << "Number of cuts: " << cuts << endl;
    cout << "Max cells possible: " << maximum_cell << endl;
   
    return 0;
}
