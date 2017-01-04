/**
 * @file   GeometryTestAlg.h
 * @brief  Unit test for geometry functionalities
 * @date   2011/02/17
 * @author brebel@fnal.gov
 * @see    GeometryTestAlg.cxx
 * 
 * Refactored by Gianluca Petrillo on May 5th, 2015.
 */

#ifndef GEO_GEOMETRYTESTALG_H
#define GEO_GEOMETRYTESTALG_H

// LArSoft includes
#include "TestUtils/NameSelector.h"

// C/C++ standard libraries
#include <string>
#include <set>
#include <vector>
#include <array>
#include <memory> // std::unique_ptr<>

// forward declarations
namespace fhicl {
  class ParameterSet;
}

namespace gar {
  namespace geo {
    
      // forward declarations
    class GeometryCore;
    
    namespace details {
      class TestTrackerClassBase;
    } // namespace details
    
    
    /** **************************************************************************
     * @brief Performs tests on the geometry as seen by Geometry service
     *
     * Configuration parameters
     * =========================
     *
     * - **ForgiveExceptions** (list of strings, default: empty): the categories
     *   of exceptions in this list are "forgiven" (non-fatal)
     * - **RunTests** (string list): marks which tests to run; each entry can be
     *   specified with a "+" or "-" prepended, to indicate to add or to remove
     *   the test from the set to be run, respectively. The directives act against
     *   a pre-existing default set that executes all tests, except the ones
     *   explicitly marked as excluded by default in the list below:
     *   + `CheckOverlaps` (not in default) perform overlap checks
     *   + `ThoroughCheck` (not in default) makes ROOT perform full geometry check
     *   + `FindVolumes`: checks it can find the volumes corresponding to world
     *   + `NearestChannel`:
     *   + `default`: represents the default set (optionally prepended by '@')
     *   + `!` (special): means to forget the tests configured so far; used as the
     *     first test name, removes the default list but leaves unchanged the
     *     default behaviour (the one specified with "+*" or "-*")
     * - **CheckForOverlaps** (boolean, default: false): equivalent to enabling
     *   `+CheckOverlaps` in `RunTests`
     */
    class GeometryTestAlg {
    public:
      explicit GeometryTestAlg(fhicl::ParameterSet const& pset);
      
        /// Virtual destructor
      virtual ~GeometryTestAlg() = default;
      
        /// Runs the test
      virtual void Setup(geo::GeometryCore const& new_geo) { geom = &new_geo; }
      
        /// Runs the test, returns a number of errors (very unlikely!)
      virtual unsigned int Run();
      
      
    private:
      gar::geo::GeometryCore const* geom;                     ///< pointer to geometry service provider
      std::set<std::string>         fNonFatalExceptions;
      
      // using as pointer just not to have to write the declaration in the header
      testing::NameSelector fRunTests; ///< test filter
      
      void printChannelSummary();
      void printVolBounds();
      void printDetDim();
      void printAllGeometry() const;
      void testFindVolumes();
      void testNearestChannel();
      void testStepping();
      
      bool shouldRunTests(std::string test_name) const;
      
      // single tests for testFindVolumes
      unsigned int testFindWorldVolumes();
      unsigned int testFindTPCvolumePaths();
      
      
    };
    
    
    namespace details {
        /// Class telling whether a test needs to be run
      class TestTrackerClassBase {
      public:
        using TestList_t = std::set<std::string>;
        
        virtual ~TestTrackerClassBase() = default;
        
          /// Returns whether the specified test should run
        virtual bool ShouldRun(std::string test_name) const = 0;
        
          /// Checks the test and records the request
        bool operator() (std::string test_name);
        
          /// Allow the specified test to run
        virtual void PleaseRunAlso(std::string test_name) = 0;
        
          /// Returns the tests that have been run
        TestList_t const& RunTests() const { return run; }
        
          /// Returns the tests that have been skipped
        TestList_t const& SkippedTests() const { return skipped; }
        
          /// Returns the tests that have been queried
        TestList_t QueriedTests() const;
        
          /// Checks that the validity of the configuration (after the fact)
        virtual bool CheckQueriesRegistry() const;
        
          /// Prints information about the configuration of the filter
        virtual void PrintConfiguration(std::ostream&) const;
        
      protected:
        TestList_t run;     ///< requested tests that should be run
        TestList_t skipped; ///< requested tests that should be skipped
        
        virtual void RecordRequest(std::string test_name, bool bRun);
        
          /// Checks the test and records the request
        virtual bool Query(std::string test_name);
        
          /// Adds a vector of tests into a test set
        static void CopyList(TestList_t                    & dest,
                             std::vector<std::string> const& from);
        
      }; // class TestTrackerClassBase
      
    } // namespace details
    
    
  } // namespace geo
} // namespace gar

#endif // GEO_GEOMETRYTESTALG_H
