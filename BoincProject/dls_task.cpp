// This file is part of BOINC.
// http://boinc.berkeley.edu
// Copyright (C) 2008 University of California
//
// BOINC is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// BOINC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with BOINC.  If not, see <http://www.gnu.org/licenses/>.

// This program serves as both
// - An example BOINC application, illustrating the use of the BOINC API
// - A program for testing various features of BOINC
//
// NOTE: this file exists as both
// boinc/apps/upper_case.cpp
// and
// boinc_samples/example_app/uc2.cpp
// If you update one, please update the other!

// The program converts a mixed-case file to upper case:
// read "in", convert to upper case, write to "out"
//
// command line options
// --cpu_time N: use about N CPU seconds after copying files
// --critical_section: run most of the time in a critical section
// --early_exit: exit(10) after 30 chars
// --early_crash: crash after 30 chars
// --run_slow: sleep 1 second after each character
// --trickle_up: sent a trickle-up message
// --trickle_down: receive a trickle-up message
// --network_usage: tell the client we used some network
//

#ifdef _WIN32
#include "boinc_win.h"
#else
#include "config.h"
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#endif

#include "str_util.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"
#include "graphics2.h"
#include "uc2.h"

#include "DiagonalLatinSquares/ComputationTask.hpp"
#include "DiagonalLatinSquares/DLS.hpp"
#include <fstream>

using std::string;

#define CHECKPOINT_FILE "upper_case_state"
#define INPUT_FILENAME "in"
#define OUTPUT_FILENAME "out"

bool run_slow = false;
bool early_exit = false;
bool early_crash = false;
bool early_sleep = false;
bool trickle_up = false;
bool trickle_down = false;
bool critical_section = false;    // run most of the time in a critical section
bool report_fraction_done = true;
bool network_usage = false;
double cpu_time = 20, comp_result;


int main(int argc, char **argv) {
    int i, retval;
    char input_path[512], output_path[512], chkpt_path[512], buf[256];

    for (i=0; i<argc; i++) {
        if (strstr(argv[i], "early_exit")) early_exit = true;
        if (strstr(argv[i], "early_crash")) early_crash = true;
        if (strstr(argv[i], "early_sleep")) early_sleep = true;
        if (strstr(argv[i], "run_slow")) run_slow = true;
        if (strstr(argv[i], "critical_section")) critical_section = true;
        if (strstr(argv[i], "network_usage")) network_usage = true;
        if (strstr(argv[i], "cpu_time")) {
            cpu_time = atof(argv[++i]);
        }
        if (strstr(argv[i], "trickle_up")) trickle_up = true;
        if (strstr(argv[i], "trickle_down")) trickle_down = true;
    }
    retval = boinc_init();
    if (retval) {
        fprintf(stderr, "%s boinc_init returned %d\n",
                boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit(retval);
    }

    fprintf(stderr, "%s app started; CPU time %f, flags:%s%s%s%s%s%s%s\n",
            boinc_msg_prefix(buf, sizeof(buf)),
            cpu_time,
            early_exit?" early_exit":"",
            early_crash?" early_crash":"",
            early_sleep?" early_sleep":"",
            run_slow?" run_slow":"",
            critical_section?" critical_section":"",
            trickle_up?" trickle_up":"",
            trickle_down?" trickle_down":""
    );

    boinc_resolve_filename(INPUT_FILENAME, input_path, sizeof(input_path));
    boinc_resolve_filename(OUTPUT_FILENAME, output_path, sizeof(output_path));
    boinc_resolve_filename(CHECKPOINT_FILE, chkpt_path, sizeof(chkpt_path));

    if (network_usage) {
        boinc_network_usage(5., 17.);
    }

    // Create computation task
    ComputationTask task(input_path, output_path);
    task.GoToCheckpoint(chkpt_path);

    while (task.DoIteration()) {
        if (boinc_time_to_checkpoint()) {
            retval = task.MakeCheckpoint(chkpt_path);
            if (retval) {
                fprintf(stderr, "%s APP: DLS checkpoint failed %d\n",
                        boinc_msg_prefix(buf, sizeof(buf)), retval
                );
                exit(retval);
            }
            boinc_checkpoint_completed();
        }
        boinc_fraction_done(task.GetFractionDone());
    }
    task.Finish();

    boinc_fraction_done(1.0);
    boinc_finish(0);

    return 0;
}
