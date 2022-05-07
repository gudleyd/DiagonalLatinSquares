// This file is part of BOINC.
// http://boinc.berkeley.edu
// Copyright (C) 2015 University of California
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

// A sample assimilator that:
// 1) if success, copy the output file(s) to a directory
// 2) if failure, append a message to an error log

#include <vector>
#include <string>
#include <cstdlib>
#include <sqlite3.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include "boinc_db.h"
#include "error_numbers.h"
#include "filesys.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "sched_config.h"
#include "assimilate_handler.h"

using std::vector;
using std::string;

const char* outdir = "results";
const char* db_path = "results.db";

struct DBCloser {
    sqlite3 **db;

    DBCloser(sqlite3 **db_): db(db_) { }

    ~DBCloser() {
        sqlite3_close(*db);
    }

};


sqlite3 *db;
char *zErrMsg = 0;
int retval;
char buf[1024];
DBCloser closer(&db);

int write_error(const char* p) {
    static FILE* f = 0;
    if (!f) {
        char path[1024];
        sprintf(path, "%s/errors", outdir);
        f = fopen(config.project_path(path), "a");
        if (!f) return ERR_FOPEN;
    }
    fprintf(f, "%s", p);
    fflush(f);
    return 0;
}

int assimilate_handler_init(int argc, char** argv) {
    retval = sqlite3_open(config.project_path(db_path), &db);
    if (retval) {
        sprintf(buf, "Cannot open db");
        return write_error(buf);
    }

    std::string sql = "CREATE TABLE IF NOT EXISTS RESULT("  \
      "ID             TEXT    NOT NULL," \
      "SEQNO          INTEGER NOT NULL," \
      "POSINSEQ       INTEGER NOT NULL," \
      "SUBRECTANGLES  INTEGER NOT NULL," \
      "FULL_CYCLES    INTEGER NOT NULL," \
      "PARTIAL_CYCLES INTEGER NOT NULL );";

    retval = sqlite3_exec(db, sql.c_str(), 0, 0, &zErrMsg);
    if (retval != SQLITE_OK) {
        sprintf(buf, "DB Table Creation Error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return write_error(buf);
    }

    for (int i=1; i<argc; i++) {
        if (!strcmp(argv[i], "--outdir")) {
            outdir = argv[++i];
        } else {
            fprintf(stderr, "bad arg %s\n", argv[i]);
        }
    }
    return 0;
}

void assimilate_handler_usage() {
    // describe the project specific arguments here
    fprintf(stderr,
            "    Custom options:\n"
            "    [--outdir X]  output dir for result files\n"
    );
}

int assimilate_handler(
        WORKUNIT& wu, vector<RESULT>& /*results*/, RESULT& canonical_result
) {
    unsigned int i;

    retval = boinc_mkdir(config.project_path(outdir));
    if (retval) return retval;

    if (wu.canonical_resultid) {
        vector<OUTPUT_FILE_INFO> output_files;
        get_output_file_infos(canonical_result, output_files);
        unsigned int n = output_files.size();
        for (i=0; i<n; i++) {
            OUTPUT_FILE_INFO& fi = output_files[i];
            std::ifstream result(fi.path.c_str(), std::ifstream::in);

            // split string to vector of tokens by _
            std::vector<std::string> tokens;
            std::stringstream ss(fi.path);
            std::string token;
            while (std::getline(ss, token, '_')) {
                tokens.push_back(token);
            }

            std::string srName, seqno = tokens.size() > 3 ?  tokens[3] : "-1";
            uint64_t subrectangles, full_cycles, partial_cycles;
            std::string currentTransaction;
            int currentPos = 0;
            while (result >> srName >> subrectangles >> full_cycles >> partial_cycles) {
                currentTransaction.append("INSERT OR IGNORE INTO RESULT(ID, SEQNO, POSINSEQ, SUBRECTANGLES, FULL_CYCLES, PARTIAL_CYCLES) VALUES ('" + srName + "', " + seqno + ", " + std::to_string(currentPos) + ", " + std::to_string(subrectangles) + ", " + std::to_string(full_cycles) + ", " + std::to_string(partial_cycles) + ");");
                currentPos += 1;
            }
            if (!currentTransaction.empty()) {
                retval = sqlite3_exec(db, currentTransaction.c_str(), 0, 0, &zErrMsg);
                if (retval != SQLITE_OK) {
                    sprintf(buf, "DB Transaction Error: %s\n", zErrMsg);
                    sqlite3_free(zErrMsg);
                    return write_error(buf);
                }
            }
        }
    } else {
        sprintf(buf, "%s: 0x%x\n", wu.name, wu.error_mask);
        return write_error(buf);
    }
    return 0;
}
