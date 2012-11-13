#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <sstream>
#include "../../src/dataio.hpp"
#include "../../src/character.hpp"
#include "../../src/genetree.hpp"
#include "../../src/utility.hpp"
#include "testutils.hpp"

#define DEFAULT_NREPS  1

class RunClock {
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> TimePointType;
    public:
        RunClock(const std::string& operation)
            : operation_(operation)
            , mean_elapsed_microseconds_(-1)
            , mean_elapsed_seconds_(-1) {
        }
        void start() {
            this->mean_elapsed_seconds_ = -1;
            this->begin_ = std::move(std::chrono::system_clock::now());
        }
        void stop() {
            TimePointType end = std::chrono::system_clock::now();
            unsigned long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - this->begin_).count();
            this->elapsed_microseconds_.push_back(ms);
        }
        double get_mean_elapsed_microseconds() {
            if (this->mean_elapsed_microseconds_ < 0) {
                unsigned long total_microseconds = 0;
                for (auto ms : this->elapsed_microseconds_) {
                    total_microseconds += ms;
                }
                this->mean_elapsed_microseconds_ = static_cast<double>(total_microseconds) / (static_cast<double>(this->elapsed_microseconds_.size()));
            }
            return this->mean_elapsed_microseconds_;
        }
        double get_mean_elapsed_seconds() {
            if (this->mean_elapsed_seconds_ < 0) {
                this->mean_elapsed_seconds_ = static_cast<double>(this->get_mean_elapsed_microseconds()) / (1e6);
            }
            return this->mean_elapsed_seconds_;
        }
        void print(std::ostream& out) {
            out << std::setw(25) << std::left << this->operation_ << " ";
            out << std::setprecision(12) << this->get_mean_elapsed_microseconds() << std::endl;
        }

    private:
        std::string                 operation_;
        TimePointType               begin_;
        std::vector<unsigned long>  elapsed_microseconds_;
        double                      mean_elapsed_microseconds_;
        double                      mean_elapsed_seconds_;
}; // RunClock

bool cmp_results(RunClock * a, RunClock * b) {
    if (a->get_mean_elapsed_seconds() < b->get_mean_elapsed_seconds()) {
        return true;
    } else {
        return false;
    }
}

class TimeLogger {

    public:
        ~TimeLogger() {
            for (auto log : this->logs_) {
                delete log;
            }
        }
        RunClock * get_timer(const std::string& operation) {
            RunClock * rc = nullptr;
            std::map<std::string, RunClock *>::iterator rci = this->run_clocks_.find(operation);
            if (rci != this->run_clocks_.end()) {
                rc = rci->second;
            } else {
                rc = new RunClock(operation);
                this->run_clocks_[operation] = rc;
                this->logs_.push_back(rc);
                this->logs_by_operation_[operation].push_back(rc);
            }
            return rc;
        }
        void summarize(std::ostream& out) {
            std::sort(this->logs_.begin(), this->logs_.end(), &cmp_results);
            for (auto rc : this->logs_) {
                rc->print(out);
            }
        }
        void summarize_by_operation(std::ostream& out) {
            for (auto log_by_operation : this->logs_by_operation_) {
                std::string operation = log_by_operation.first;
                out << "\n\n### " << operation << " ###\n\n";
                std::vector<RunClock *> logs = log_by_operation.second;
                std::sort(logs.begin(), logs.end(), &cmp_results);
                for (auto rc : logs) {
                    rc->print(out);
                }
            }
        }
        void summarize_best_by_operation(std::ostream& out) {
            for (auto log_by_operation : this->logs_by_operation_) {
                std::string operation = log_by_operation.first;
                std::vector<RunClock *> logs = log_by_operation.second;
                std::sort(logs.begin(), logs.end(), &cmp_results);
                (*logs.begin())->print(out);
            }
        }

    private:
        std::vector<RunClock *>                        logs_;
        std::map<std::string, std::vector<RunClock *>> logs_by_operation_;
        std::map<std::string, RunClock *>              run_clocks_;

}; // TimeLogger

void run_tree_construction(const std::string& tree_filepath,
        TimeLogger& time_logger,
        unsigned long nreps) {
    RunClock * clock = nullptr;
    for (unsigned long i = 0; i < nreps; ++i) {
        std::vector<treeshrew::GeneTree *> trees;
        clock = time_logger.get_timer("Tree Construction");
        clock->start();
        treeshrew::treeio::read_from_filepath(trees, tree_filepath, "newick");
        clock->stop();
    }
}

void run_postorder_iteration(const std::vector<treeshrew::GeneTree *>& trees,
        TimeLogger& time_logger,
        unsigned long nreps) {
    RunClock * clock = nullptr;
    for (unsigned long i = 0; i < nreps; ++i) {
        for (auto & tree : trees) {
            clock = time_logger.get_timer("Postorder Iteration");
            clock->start();
            for (treeshrew::GeneTreeNode::postorder_iterator ndi = tree->postorder_begin(); ndi != tree->postorder_end(); ++ndi) {
            }
            clock->stop();
        }
    }
}

void run_leaf_iteration(const std::vector<treeshrew::GeneTree *>& trees,
        TimeLogger& time_logger,
        unsigned long nreps) {
    RunClock * clock = nullptr;
    for (unsigned long i = 0; i < nreps; ++i) {
        for (auto & tree : trees) {
            clock = time_logger.get_timer("Leaf Iteration");
            clock->start();
            for (treeshrew::GeneTreeNode::leaf_iterator ndi = tree->leaf_begin(); ndi != tree->leaf_end(); ++ndi) {
            }
            clock->stop();
        }
    }
}

void run_child_iteration(const std::vector<treeshrew::GeneTree *>& trees,
        TimeLogger& time_logger,
        unsigned long nreps) {
    RunClock * clock = nullptr;
    for (auto & tree : trees) {
        clock = time_logger.get_timer("Leaf Iteration");
        clock->start();
        for (treeshrew::GeneTreeNode::postorder_iterator ndi = tree->postorder_begin(); ndi != tree->postorder_end(); ++ndi) {
            auto node = *ndi;
            for (treeshrew::GeneTreeNode::child_iterator chi = node->children_begin();
                    chi != node->children_end();
                    ++chi) {
            }
        }
        clock->stop();
    }
}

int main(int argc, char * argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <NEWICK-TREEFILE> <FASTA-DATAFILE>" << std::endl;
        exit(1);
    }
    std::string tree_filepath(argv[1]);
    std::string data_filepath(argv[2]);
    unsigned long nreps = DEFAULT_NREPS;
    if (argc >= 4) {
        std::istringstream(argv[3]) >> nreps;
    }

    treeshrew::NucleotideSequences data;
    treeshrew::sequenceio::read_from_filepath(data, data_filepath, "fasta");
    std::vector<treeshrew::GeneTree *> trees;
    treeshrew::treeio::read_from_filepath(trees, tree_filepath, "newick");

    TimeLogger time_logger;
    run_tree_construction(tree_filepath, time_logger, nreps);
    run_postorder_iteration(trees, time_logger, nreps);
    run_leaf_iteration(trees, time_logger, nreps);
    run_child_iteration(trees, time_logger, nreps);

    time_logger.summarize(std::cout);
}


