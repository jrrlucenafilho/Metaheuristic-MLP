#include "Data.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

typedef struct {
    vector<int> sequence;
    double cost = 0.0;
} Solution;

//Aux structures. Holds relevant data for a given subsequence MLP
struct Subsequence {
    double duration, accumulatedCost;   //T, C
    int delayCost; //W
    int firstNode, lastNode; //pertaining to subsequence

    inline static Subsequence Concatenate(Subsequence& sigma1, Subsequence& sigma2, double** m){
        Subsequence sigma;
        double temp = m[sigma1.lastNode][sigma2.firstNode];

        sigma.delayCost = sigma1.delayCost + sigma2.delayCost;
        sigma.duration = sigma1.duration + temp + sigma2.duration;
        sigma.accumulatedCost = sigma1.accumulatedCost + sigma2.delayCost * (sigma1.duration + temp) + sigma2.accumulatedCost;
    
        sigma.firstNode = sigma1.firstNode;
        sigma.lastNode = sigma2.lastNode;

        return sigma;
    }
};

typedef struct {
    clock_t accumulatedTime = 0;
    clock_t beginTime = 0;
    clock_t endTime = 0;
} my_time_t;

//Timers, they'll hold the sum of time in clock cycles over all iterations (for avg calc)
my_time_t buildSol_time, bestImprovSwap_time, orOpt_time, orOpt2_time, orOpt3_time, two_opt_time, disturbance_time;

my_time_t* buildSol_time_ptr = &buildSol_time; 
my_time_t* bestImprovSwap_time_ptr = &bestImprovSwap_time;
my_time_t* orOpt_time_ptr = &orOpt_time;
my_time_t* orOpt2_time_ptr = &orOpt2_time;
my_time_t* orOpt3_time_ptr = &orOpt3_time;
my_time_t* two_opt_time_ptr = &two_opt_time;
my_time_t* disturbance_time_ptr = &disturbance_time;


//Prints the time for each Neighborhood Structure
void PrintNBSTimers()
{
    cout << "\tBuildSolution time: " << buildSol_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tBestImprovementSwap time: " << bestImprovSwap_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tOrOpt time: " << orOpt_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tOrOpt2 time: " << orOpt2_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tOrOpt3 time: " << orOpt3_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\t2-Opt time: " << two_opt_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
    cout << "\tDisturbance time: " << disturbance_time_ptr->accumulatedTime / static_cast<double>(CLOCKS_PER_SEC) << " s\n";
}

/**
 * @brief Inits all subseq values of a given param solution. Where ele seq_matrix[i][j] holds all subseq info in solution seq
 * @var solNodeQuant Number of nodes in current instance
 * @param solution
 * @param subseqMatrix Matrix holding all subseq info in solution seq
**/
void UpdateAllSubSeqs(Solution& solution, vector<vector<Subsequence>>& subseqMatrix, double** m)
{
    int solNodeQuant = solution.sequence.size();

    //First, deals with subseqs comprised of only 1 node
    for(int i = 0; i < solNodeQuant; i++){
        subseqMatrix[i][i].delayCost = (i > 0);
        subseqMatrix[i][i].accumulatedCost = 0;
        subseqMatrix[i][i].duration = 0;
        subseqMatrix[i][i].firstNode = solution.sequence[i];
        subseqMatrix[i][i].lastNode = solution.sequence[i];
    }

    //Now deal with the non-single-node subseqs
    for(int i = 0; i < solNodeQuant; i++){
        for(int j = i + 1; j < solNodeQuant; j++){
            subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
        }
    }

    //Now deals with inverted subseqs (for 2-opt)
    for(int i = solNodeQuant - 1; i >= 0; i--){
        for(int j = i - 1; j >= 0; j--){
            subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j + 1], subseqMatrix[j][j], m);
        }
    }
}

double CalculateSequenceCost(vector<int>& sequence, double** m, int dimension)
{
    double sum = 0;
    double branchTotalCost = 0;   //Holds the accumulated branch cost (including current node)
    double formerNodeCost = 0;

    for(int i = 0, j = 1; i < dimension; i++, j++){
        formerNodeCost = m[sequence[i]][sequence[j]];   //Cost of the edge between the former and current node
        sum += formerNodeCost;
        branchTotalCost += sum;    //holds all costs of previous edges on the branch
    }
    
    return branchTotalCost;
}

void RecalculateCL(vector<pair<double, int>>& candidatesList, double** m, int dimension, int selectedNode)
{
    for(int i = 1; i <= dimension; i++){
        //First checks if the node is in the list
        auto iter = find_if(candidatesList.begin(), candidatesList.end(), [i](const pair<double, int>& element)
                          { return element.second == i; });

        //Skips if node is itself or if it isn't in CL
        if((i == selectedNode) || iter == candidatesList.end()){
            continue;
        }

        //If it is in CL, gets the cost from the matrix, in ref to current selectedNode node
        iter->first = m[i][selectedNode];
    }
}

/**
 * @brief Builds a fair solution (though still far from optimized) Using Greedy Randomized Adaptive Search (Insertion-by-cheapest)
 * @return Solution
 */
Solution BuildSolution(double** distMatrix, int dimension)
{
    buildSol_time_ptr->beginTime = std::clock();

    Solution solution;
    vector<pair<double, int>> candidatesList;   //pair <cost, node>
    int selected = -1;  //Just so it doesn't recalc first iter (doesn't need to)

    solution.sequence.reserve(dimension + 1);

    //First begin building seq (already inits at cost being zero)
    solution.sequence = {1};

    //Populating CL, calc'ing in ref to first node's distance to each
    for(int i = 2; i <= dimension; i++){
        candidatesList.push_back({distMatrix[1][i], i});
    }

    while(!candidatesList.empty()){
        //Changing ref point to current, recalc'ing costs regarding last inserted node (if not first iter, which begins on first node (row 1))
        if(selected != -1){
            RecalculateCL(candidatesList, distMatrix, dimension, candidatesList[selected].second);
        }

        //Sorts the list of candidates by cost
        sort(candidatesList.begin(), candidatesList.end());

        //Somewhat random
        double alpha = ((double)rand() + 1) / RAND_MAX;
        selected = rand() % ((int)ceil(alpha * candidatesList.size()));

        //Add it to initial sol seq and remove from CL
        solution.sequence.push_back(candidatesList[selected].second);
        candidatesList.erase(candidatesList.begin() + selected);
    }

    solution.sequence.push_back(1);

    buildSol_time_ptr->endTime = std::clock();
    buildSol_time_ptr->accumulatedTime += buildSol_time_ptr->endTime - buildSol_time_ptr->beginTime;

    return solution;
}

/*BestImprovement funcs with differing methods to be used for RVND in LocalSearch()*/
//Tries to get the best swap possible in the solution
//As in the one that best minimizes the solution's cost
//"m" here means distMatrix
bool BestImprovementSwap(Solution& solution, vector<vector<Subsequence>>& subseqMatrix, double** m, int dimension)
{
    bestImprovSwap_time_ptr->beginTime = std::clock();

    //Inits
    Subsequence sigma1, sigma2, sigma3, sigmaLast;
    double bestCost = subseqMatrix[0][dimension].accumulatedCost;
    int graphSize = dimension + 1;
    int best_i = 0, best_j = 0;

    for(int i = 1; i < graphSize - 1; i++){
        for(int j = i + 1; j < graphSize - 1; j++){
            if(j == i + 1){ //Case of neightboring nodes
                sigma1 = Subsequence::Concatenate(subseqMatrix[0][i - 1], subseqMatrix[j][j], m);
                sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[i][i], m);
                sigmaLast = Subsequence::Concatenate(sigma2, subseqMatrix[j + 1][dimension], m);
            }else{
                sigma1 = Subsequence::Concatenate(subseqMatrix[0][i - 1], subseqMatrix[j][j], m);
                sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[i + 1][j - 1], m);
                sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[i][i], m);
                sigmaLast = Subsequence::Concatenate(sigma3, subseqMatrix[j + 1][dimension], m);
            }

            if(sigmaLast.accumulatedCost < bestCost){
                bestCost = sigmaLast.accumulatedCost;
                best_i = i;
                best_j = j;
            }
        }
    }

    //Actual swap, only happens when overall solution cost lowers
    if(bestCost < solution.cost){
        solution.cost = bestCost;
        swap(solution.sequence[best_i], solution.sequence[best_j]);
        swap(subseqMatrix[best_i][best_i], subseqMatrix[best_j][best_j]);

        //Concats
        for(int i = 0; i < best_i; i++){
            for(int j = best_i; j < graphSize; j++){
                subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
            }
        }

        for(int i = best_i; i <= best_j; i++){  
            for(int j = i + 1; j < graphSize; j++){ 
                subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
            }
        }

        bestImprovSwap_time_ptr->endTime = std::clock();
        bestImprovSwap_time_ptr->accumulatedTime += bestImprovSwap_time_ptr->endTime - bestImprovSwap_time_ptr->beginTime;
        
        return true;
    }

    bestImprovSwap_time_ptr->endTime = std::clock();
    bestImprovSwap_time_ptr->accumulatedTime += bestImprovSwap_time_ptr->endTime - bestImprovSwap_time_ptr->beginTime;

    return false;
}

bool BestImprovement2Opt(Solution& solution, vector<vector<Subsequence>>& subseqMatrix, double** m, int dimension)
{
    two_opt_time_ptr->beginTime = std::clock();

    Subsequence sigma1, sigma2;
    double bestCost = subseqMatrix[0][dimension].accumulatedCost;
    int graphSize = dimension + 1;
    int best_i = 0, best_j = 0;

    for(int i = 1; i < graphSize - 1; i++){
        for(int j = i + 1; j < graphSize - 1; j++){
            sigma1 = Subsequence::Concatenate(subseqMatrix[0][i - 1], subseqMatrix[j][i], m);
            sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[j + 1][dimension], m);

            if(sigma2.accumulatedCost < bestCost){
                bestCost = sigma2.accumulatedCost;
                best_i = i;
                best_j = j;
            }
        }
    }

    //Actual swap
    if(bestCost < solution.cost){
        reverse(solution.sequence.begin() + best_i, solution.sequence.begin() + best_j + 1);

        solution.cost = bestCost;
        UpdateAllSubSeqs(solution, subseqMatrix, m);

    for(int i = best_i; i <= best_j; i++){
        subseqMatrix[i][i].delayCost = (i > 0);
        subseqMatrix[i][i].accumulatedCost = 0;
        subseqMatrix[i][i].duration = 0;
        subseqMatrix[i][i].firstNode = solution.sequence[i];
        subseqMatrix[i][i].lastNode = solution.sequence[i];
    }

        //Concats
        for(int i = best_i; i <= best_j; i++){
            for(int j = i + 1; j < graphSize; j++){
                subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
            }
        }

        for(int i = best_i; i <= best_j; i++){
            for(int j = i + 1; j < graphSize; j++){
                subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
            }
        }

        two_opt_time_ptr->endTime = std::clock();
        two_opt_time_ptr->accumulatedTime += two_opt_time_ptr->endTime - two_opt_time_ptr->beginTime;

        return true;
    }

    two_opt_time_ptr->endTime = std::clock();
    two_opt_time_ptr->accumulatedTime += two_opt_time_ptr->endTime - two_opt_time_ptr->beginTime;

    return false;
}

bool BestImprovementOrOpt(Solution& solution, vector<vector<Subsequence>>& subseqMatrix, double** m, int dimension, int movedBlockSize)
{
    //Just to account for var creation time, smallish
    my_time_t var_creation_time;
    var_creation_time.beginTime = std::clock();

    Subsequence sigma1, sigma2, sigma3;
    double bestCost = subseqMatrix[0][dimension].accumulatedCost;
    int graphSize = dimension + 1;
    int best_i = 0, best_j = 0;

    var_creation_time.endTime = std::clock();
    var_creation_time.accumulatedTime += var_creation_time.endTime - var_creation_time.beginTime;

    //Reinsertion Case
    if(movedBlockSize == 1){
        orOpt_time_ptr->beginTime = std::clock();

        for(int i = 1; i < graphSize - 1; i++){        
            for(int j = 1; j < graphSize - 1; j++){
                if(i != j){
                    if(i < j){
                        sigma1 = Subsequence::Concatenate(subseqMatrix[0][i - 1], subseqMatrix[i + 1][j], m);
                        sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[i][i], m);
                        sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[j + 1][dimension], m);
                    }else{
                        sigma1 = Subsequence::Concatenate(subseqMatrix[0][j - 1], subseqMatrix[i][i], m);
                        sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[j][i - 1], m);
                        sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[i + 1][dimension], m);
                    }

                    if(sigma3.accumulatedCost < bestCost){
                        bestCost = sigma3.accumulatedCost;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        if(bestCost < solution.cost){
            int reinsertPortion = solution.sequence[best_i];
            solution.sequence.erase(solution.sequence.begin() + best_i);
            solution.sequence.insert(solution.sequence.begin() + best_j, reinsertPortion);

            solution.cost = bestCost;

            //Swap and concats
            if(best_i > best_j){
                swap(best_i, best_j);
            }

            for(int i = best_i; i <= best_j; i++){
                subseqMatrix[i][i].delayCost = (i > 0);
                subseqMatrix[i][i].accumulatedCost = 0;
                subseqMatrix[i][i].duration = 0;
                subseqMatrix[i][i].firstNode = solution.sequence[i];
                subseqMatrix[i][i].lastNode = solution.sequence[i];
            }

            for(int i = 0; i < best_i; i++){  
                for(int j = best_i; j < graphSize; j++){
                    subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                    subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
                }
            }

            for(int i = best_i; i <= best_j; i++){  
                for(int j = i + 1; j < graphSize; j++){
                    subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m); 
                    subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);   
                }
            }

            orOpt_time_ptr->endTime = std::clock();
            orOpt_time_ptr->accumulatedTime += orOpt_time_ptr->endTime - orOpt_time_ptr->beginTime;
            orOpt_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;

            return true;
        }

        orOpt_time_ptr->endTime = std::clock();
        orOpt_time_ptr->accumulatedTime += orOpt_time_ptr->endTime - orOpt_time_ptr->beginTime;
        orOpt_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;

        return false;
    }

    //OrOpt-2 case
    if(movedBlockSize == 2){
        orOpt2_time_ptr->beginTime = std::clock();

        for(int i = 1; i < graphSize - 2; i++){
            for(int j = 1; j < graphSize - 3; j++){
                if(i != j){
                    if(i < j){
                        sigma1 = Subsequence::Concatenate(subseqMatrix[0][i - 1], subseqMatrix[i + 2][j + 1], m);
                        sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[i][i + 1], m);
                        sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[j + 2][dimension], m);
                    }else{
                        sigma1 = Subsequence::Concatenate(subseqMatrix[0][j - 1], subseqMatrix[i][i + 1], m);
                        sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[j][i - 1], m);
                        sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[i + 2][dimension], m);
                    }

                    if(sigma3.accumulatedCost < bestCost){
                        bestCost = sigma3.accumulatedCost;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        if(bestCost < solution.cost){
            vector<int> reinsertPortion(solution.sequence.begin() + best_i, solution.sequence.begin() + best_i + 2);
            solution.sequence.erase(solution.sequence.begin() + best_i, solution.sequence.begin() + best_i + 2);
            solution.sequence.insert(solution.sequence.begin() + best_j, reinsertPortion.begin(), reinsertPortion.end());

            solution.cost = bestCost;

            if(best_i > best_j){
                swap(best_i, best_j);
            }

            for(int i = best_i; i <= best_j + 1; i++){
                subseqMatrix[i][i].delayCost = (i > 0);
                subseqMatrix[i][i].accumulatedCost = 0;
                subseqMatrix[i][i].duration = 0;
                subseqMatrix[i][i].firstNode = solution.sequence[i];
                subseqMatrix[i][i].lastNode = solution.sequence[i];
            }

            for(int i = 0; i < best_i; i++){  
                for(int j = best_i; j < graphSize; j++){
                    subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                    subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
                }
            }

            for(int i = best_i; i <= best_j + 1; i++){  
                for(int j = i + 1; j < graphSize; j++){
                    subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m); 
                    subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);   
                }
            }

            orOpt2_time_ptr->endTime = std::clock();
            orOpt2_time_ptr->accumulatedTime += orOpt2_time_ptr->endTime - orOpt2_time_ptr->beginTime;
            orOpt2_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;

            return true;
        }
        
        orOpt2_time_ptr->endTime = std::clock();
        orOpt2_time_ptr->accumulatedTime += orOpt2_time_ptr->endTime - orOpt2_time_ptr->beginTime;
        orOpt_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;

        return false;
    }

    //OrOpt-3 case
    if(movedBlockSize == 3){
        orOpt3_time_ptr->beginTime = std::clock();

        for(int i = 1; i < graphSize - 3; i++){
            for(int j = 1; j < graphSize - 4; j++){
                if(i != j){
                    if(i < j){
                        sigma1 = Subsequence::Concatenate(subseqMatrix[0][i - 1], subseqMatrix[i + 3][j + 2], m);
                        sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[i][i + 2], m);
                        sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[j + 3][dimension], m);
                    }else{
                        sigma1 = Subsequence::Concatenate(subseqMatrix[0][j - 1], subseqMatrix[i][i + 2], m);
                        sigma2 = Subsequence::Concatenate(sigma1, subseqMatrix[j][i - 1], m);
                        sigma3 = Subsequence::Concatenate(sigma2, subseqMatrix[i + 3][dimension], m);
                    }

                    if(sigma3.accumulatedCost < bestCost){
                        bestCost = sigma3.accumulatedCost;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }

        //Concats
        if((bestCost < solution.cost)){
            vector<int> reinsertSequence(solution.sequence.begin() + best_i, solution.sequence.begin() + best_i + 3);
            solution.sequence.erase(solution.sequence.begin() + best_i, solution.sequence.begin() + best_i + 3);
            solution.sequence.insert(solution.sequence.begin() + best_j, reinsertSequence.begin(), reinsertSequence.end());

            solution.cost = bestCost;

            if(best_i > best_j){
                swap(best_i, best_j);
            }

            for(int i = best_i; i <= best_j + 2; i++){
                subseqMatrix[i][i].delayCost = (i > 0);
                subseqMatrix[i][i].accumulatedCost = 0;
                subseqMatrix[i][i].duration = 0;
                subseqMatrix[i][i].firstNode = solution.sequence[i];
                subseqMatrix[i][i].lastNode = solution.sequence[i];
            }

            for(int i = 0; i < best_i; i++){  
                for(int j = best_i; j < graphSize; j++){
                    subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m);
                    subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);
                }
            }

            for(int i = best_i; i <= best_j + 2; i++){  
                for(int j = i + 1; j < graphSize; j++){
                    subseqMatrix[i][j] = Subsequence::Concatenate(subseqMatrix[i][j - 1], subseqMatrix[j][j], m); 
                    subseqMatrix[j][i] = Subsequence::Concatenate(subseqMatrix[j][j], subseqMatrix[j - 1][i], m);   
                }
            }

            orOpt3_time_ptr->endTime = std::clock();
            orOpt3_time_ptr->accumulatedTime += orOpt3_time_ptr->endTime - orOpt3_time_ptr->beginTime;
            orOpt3_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;

            return true;
        }

        orOpt3_time_ptr->endTime = std::clock();
        orOpt3_time_ptr->accumulatedTime += orOpt3_time_ptr->endTime - orOpt3_time_ptr->beginTime;
        orOpt_time_ptr->accumulatedTime += var_creation_time.accumulatedTime;
    
        return false;
    }

    return false;
}

//Increases the quality of current iteration's solution
//Does that by modifying the solution and evaluing the impact of each change
//Using the Random Variable Neighborhood Descent method
//Which just tests different neighborhood structures with a tad of randomness when choosing
//discarding whichever makes cost higher than currCost
void LocalSearch(Solution& solution, vector<vector<Subsequence>>& subseqMatrix, double** distMatrix, int dimension)
{
    vector<int> NH_structures = {1, 2, 3, 4, 5};
    bool solutionImproved = false;

    while(!NH_structures.empty()){
        int rand_n = rand() % NH_structures.size();

        switch(NH_structures[rand_n]){
            case 1:
                solutionImproved = BestImprovementSwap(solution, subseqMatrix, distMatrix, dimension);
                break;
            case 2: 
                solutionImproved = BestImprovement2Opt(solution, subseqMatrix, distMatrix, dimension);
                break;
            case 3:
                solutionImproved = BestImprovementOrOpt(solution, subseqMatrix, distMatrix, dimension, 1);
                break;
            case 4:
                solutionImproved = BestImprovementOrOpt(solution, subseqMatrix, distMatrix, dimension, 2);
                break;
            case 5:
                solutionImproved = BestImprovementOrOpt(solution, subseqMatrix, distMatrix, dimension, 3);
                break;
        }

        //Checks if solution's improved, erasing neighborhood structure from vec if not
        if(solutionImproved){
            NH_structures = {1, 2, 3, 4, 5};
        }else{
            NH_structures.erase(NH_structures.begin() + rand_n);
        }
    }
}

int BoundedRand(int min, int max)
{
    return min + rand() % (max - min + 1);
}

//Return Solution, may have to call a cost-calc function inside here
Solution Disturbance(Solution& solution, double** m, int dimension)
{
    disturbance_time_ptr->beginTime = std::clock();

    vector<int> copiedSeq = solution.sequence;
    int subseqMaxLength = ceil(dimension / 10.0);
    Solution disturbedSol;

    //Will mark the index of first and last elements of each subsequence
    //Used to make it so they don't overlap
    int subseq1Index_begin, subseq2Index_begin;
    int subseq1Index_end, subseq2Index_end;

    //Length of each subsequence and of the space among them
    int subseq2Length, inbetweenSubseqsLength; //subseq2Length uneeded

    subseq1Index_begin = BoundedRand(1, dimension - subseqMaxLength - subseqMaxLength - 1);
    subseq1Index_end = BoundedRand(subseq1Index_begin + 1, subseq1Index_begin + subseqMaxLength - 1);
    subseq2Index_begin = BoundedRand(subseq1Index_end + 1, dimension - subseqMaxLength);
    subseq2Index_end = BoundedRand(subseq2Index_begin + 1, subseq2Index_begin + subseqMaxLength - 1);

    //Actually making the subsequences from rand calc'd indexes
    vector<int> subseq1(copiedSeq.begin() + subseq1Index_begin, copiedSeq.begin() + subseq1Index_end);
    vector<int> subseq2(copiedSeq.begin() + subseq2Index_begin, copiedSeq.begin() + subseq2Index_end);

    //Lengths and space calc
    subseq2Length = subseq2Index_end - subseq2Index_begin;
    inbetweenSubseqsLength = subseq2Index_begin - subseq1Index_end;

    //Disturbing copied sequence
    //Swapping these not-necessarily-equally-sized subseqs in the main sequence
    copiedSeq.erase(copiedSeq.begin() + subseq1Index_begin, copiedSeq.begin() + subseq1Index_end);
    copiedSeq.insert(copiedSeq.begin() + subseq1Index_begin, subseq2.begin(), subseq2.end());

    copiedSeq.erase(copiedSeq.begin() + subseq1Index_begin + subseq2Length + inbetweenSubseqsLength, copiedSeq.begin() + subseq1Index_begin + subseq2Length + inbetweenSubseqsLength + subseq2Length);
    copiedSeq.insert(copiedSeq.begin() + subseq1Index_begin + subseq2Length + inbetweenSubseqsLength, subseq1.begin(), subseq1.end());

    //Calculating cost and attributing to disturbed solution
    disturbedSol.cost = CalculateSequenceCost(copiedSeq, m, dimension);
    disturbedSol.sequence = copiedSeq;

    disturbance_time_ptr->endTime = std::clock();
    disturbance_time_ptr->accumulatedTime += disturbance_time_ptr->endTime - disturbance_time_ptr->beginTime;

    return disturbedSol;
}

Solution IteratedLocalSearch(int maxIters, int maxIterILS, Data& data)
{
    Solution bestOfAllSolution, currIterSolution, currBestSolution;
    bestOfAllSolution.cost = INFINITY;

    for(int i = 0; i < maxIters; i++){
        //Builds a beginning solution based on fair guesses
        currIterSolution = BuildSolution(data.getMatrixCost(), data.getDimension());

        vector<vector<Subsequence>> subseqMatrix(data.getDimension() + 1, vector<Subsequence>(data.getDimension() + 1));

        UpdateAllSubSeqs(currIterSolution, subseqMatrix, data.getMatrixCost());
        currBestSolution.sequence = currIterSolution.sequence;

        currIterSolution.cost = subseqMatrix[0][data.getDimension()].accumulatedCost;
        currBestSolution.cost = currIterSolution.cost;

        int iterILS = 0;

        while(iterILS <= maxIterILS){
            //Tries to enhance the fairly-guessed solution
            //By doing small modifications to it
            LocalSearch(currIterSolution, subseqMatrix, data.getMatrixCost(), data.getDimension());

            if(currIterSolution.cost < currBestSolution.cost){
                currBestSolution.cost = currIterSolution.cost;
                currBestSolution.sequence = currIterSolution.sequence;
                iterILS = 0;
            }

            //If not possible to make it better, shake the current best solution up a lil'
            //to see if we didn't just go into a 'local best pitfall'
            currIterSolution = Disturbance(currBestSolution, data.getMatrixCost(), data.getDimension());

            //Update subseqs and check if cost got better after doing disturbance
            UpdateAllSubSeqs(currIterSolution, subseqMatrix, data.getMatrixCost());
            currIterSolution.cost = CalculateSequenceCost(currIterSolution.sequence, data.getMatrixCost(), data.getDimension());

            iterILS++;
        }

        if(currBestSolution.cost < bestOfAllSolution.cost){
            bestOfAllSolution.sequence = currBestSolution.sequence;
            bestOfAllSolution.cost = currBestSolution.cost;
        }
    }

    return bestOfAllSolution;
}

int main(int argc, char** argv)
{
    int maxIter = 10;
    int maxIterILS;
    double costsSum = 0;
    auto data = Data(argc, argv[1]);
    Solution solution, finalSol;
    my_time_t ILS_iter_time;

    data.read();
    data.reformatMatrix();
    //data.printMatrixDist();
    size_t n = data.getDimension();

    cout << "Dimension: " << n << '\n';
    cout << "Wait for it...\n";

    maxIterILS = min(100, data.getDimension());

    cout << "-------------------------------\n";
    cout << "10 Iteration costs:\n";

    //10 execs for the sum avg calc
    for(int i = 0; i < 10; i++){
        srand(static_cast<unsigned int>(time(0)));

        ILS_iter_time.beginTime = std::clock();

        solution = IteratedLocalSearch(maxIter, maxIterILS, data);

        solution.cost = CalculateSequenceCost(solution.sequence, data.getMatrixCost(), data.getDimension());
        costsSum += solution.cost;

        ILS_iter_time.endTime = std::clock();
        ILS_iter_time.accumulatedTime += ILS_iter_time.endTime - ILS_iter_time.beginTime;

        cout << "Cost of s (iter " << i + 1 << "): " << solution.cost << '\n';
    }

    //Final avgs and result
    cout << "-------------------------------\n";
    cout << "Instance Name: " << data.getInstanceName() << '\n';
    cout << "Solution s = ";
    for(size_t i = 0; i < solution.sequence.size() - 1; i++){
        cout << solution.sequence[i] << " -> ";
    }
    cout << "1\n";

    cout << "Average s cost: " << costsSum / 10 << '\n';
    cout << "Average CPU execution time: " << (ILS_iter_time.accumulatedTime / static_cast<double>(CLOCKS_PER_SEC)) / 10 << " s \n";
    cout << "-------------------------------\n";
    cout << "Average Times on each Neighborhood Structure:\n";
    PrintNBSTimers();

    return 0;
}
