# Implementation Strategy: Building DAAadvisor with Claude

## 1. Project Setup & Organization

### A. Create Your Project Structure
```bash
# Create the package skeleton
mkdir DAAadvisor
cd DAAadvisor

# Create directory structure
mkdir -p R man tests/testthat inst/extdata vignettes data src

# Initialize git
git init
git add .
git commit -m "Initial package structure"
```

### B. Organize Claude Artifacts into Files
```
DAAadvisor/
├── R/
│   ├── data_profiler.R         # From artifact: daa-advisor-r-package
│   ├── method_selector.R        # From artifact: daa-method-selection-algorithm
│   ├── method_wrappers.R        # From artifact: daa-method-wrapper-template
│   ├── comparison_metrics.R     # From artifact: daa-comparison-framework
│   ├── visualization.R          # From artifact: daa-visualization-module
│   └── utils.R                  # Helper functions
├── DESCRIPTION                  # From artifact: daa-advisor-r-package
└── README.md                    # From artifact: daa-requirements-readme
```

## 2. Iterative Development Workflow with Claude

### Phase 1: Core Infrastructure (Week 1-2)
```markdown
# Prompt Template for Claude:
"I'm implementing the data_profiler.R module for DAAadvisor. 
Here's what I have so far: [paste current code]

I need to implement the following helper functions:
- .calculateSparsity()
- .assessZeroInflation()
- .detectCompositionalBias()

The functions should follow this pattern from our design: [paste relevant design section]

Please provide complete implementations with:
1. Proper error handling
2. Input validation
3. Unit tests
4. Documentation"
```

### Phase 2: Method Wrappers (Week 3-4)
```markdown
# For each method wrapper:
"I need to implement the [MethodName]Wrapper class for DAAadvisor.

Requirements:
- Inherits from DAAMethodWrapper base class
- Implements run(), check_requirements(), and validate_input()
- Returns standardized output format
- Handles edge cases (empty results, convergence failures)

Here's the base class: [paste code]
Here's an example wrapper (ALDEx2): [paste code]

Please implement the wrapper with full error handling and documentation."
```

### Phase 3: Testing Each Component
```markdown
# Testing prompt:
"I've implemented [component name]. Please help me:
1. Write comprehensive unit tests using testthat
2. Create edge case tests
3. Add integration tests with example data
4. Suggest improvements based on test coverage

Here's my implementation: [paste code]"
```

## 3. Best Practices for Claude Sessions

### A. Context Management
```r
# Start each session with context
# session_context.R
library(DAAadvisor)

# Current focus
current_module <- "method_selector"
completed_modules <- c("data_profiler", "utils")
next_steps <- c("implement scoring functions", "add ML integration")

# Key design decisions
design_principles <- list(
  sparsity_threshold = 0.7,
  min_sample_size = 5,
  confidence_calculation = "multi-factor"
)

# Paste this at the start of each Claude session
```

### B. Incremental Building
```markdown
# Good prompt structure:
1. Context: "I'm working on DAAadvisor's method selection module"
2. Previous work: "I've completed X, Y, Z functions"
3. Current goal: "Now I need to implement the scoring algorithm"
4. Specific requirements: "It should handle these edge cases..."
5. Code style: "Follow tidyverse style guide, use roxygen2 docs"
```

### C. Code Review Pattern
```markdown
# After implementing a module:
"Please review this implementation of [module]:
[paste code]

Check for:
1. Logical errors or edge cases
2. Performance optimizations
3. Better R idioms
4. Missing error handling
5. Documentation completeness"
```

## 4. Development Roadmap

### Week 1-2: Foundation
```r
# Focus on core data structures
- [ ] Implement DataProfile S3 class
- [ ] Complete all sparsity metrics
- [ ] Add compositional bias detection
- [ ] Create basic visualization methods
- [ ] Write unit tests for profiler

# Claude prompts:
"Help me implement the S3 class for DataProfile with print, summary, and plot methods"
"Create comprehensive sparsity metrics following the design in our framework"
```

### Week 3-4: Method Integration
```r
# Implement method wrappers
- [ ] Base class DAAMethodWrapper
- [ ] ALDEx2 wrapper (complete implementation)
- [ ] ANCOM-BC wrapper
- [ ] DESeq2 wrapper with zero-inflation handling
- [ ] Standardized output format

# Claude prompts:
"Implement ANCOM-BC wrapper following our ALDEx2 example"
"Add batch effect handling to DESeq2 wrapper"
```

### Week 5-6: Intelligence Layer
```r
# Method selection algorithm
- [ ] Rule-based scoring system
- [ ] Historical performance database structure
- [ ] Confidence calculations
- [ ] Explanation generation

# Claude prompts:
"Implement the rule-based scoring system from our algorithm design"
"Create explanation templates that are user-friendly"
```

### Week 7-8: Comparison Framework
```r
# Results comparison
- [ ] Consensus methods (rank product, Borda count)
- [ ] Stability assessment
- [ ] Biological validation
- [ ] Interactive visualizations

# Claude prompts:
"Implement rank product consensus following the mathematical definition"
"Create interactive plotly visualizations for method comparison"
```

## 5. Specific Claude Usage Tips

### A. For Complex Algorithms
```markdown
# Break down complex functions:
"I need to implement .calculatePhylogeneticSignal(). 
Let's break this into steps:
1. First, help me validate the phylogeny input
2. Then, calculate Pagel's lambda
3. Finally, add permutation testing

Here's step 1 context: [details]"
```

### B. For R-Specific Issues
```markdown
# R package development questions:
"I'm getting this NOTE in R CMD check: [error message]
My NAMESPACE file has: [content]
My function uses: [dependencies]
How do I properly declare these dependencies?"
```

### C. For Integration Challenges
```markdown
# When integrating multiple modules:
"I have these three modules working independently:
- data_profiler.R: [brief description]
- method_selector.R: [brief description]
- method_wrappers.R: [brief description]

Help me create the main workflow function that coordinates them"
```

## 6. Testing Strategy with Claude

### A. Test-Driven Development
```r
# First ask for tests:
"Before implementing .detectCompositionalBias(), 
help me write comprehensive tests that cover:
- Normal cases
- Edge cases (all zeros, single sample)
- Expected outputs
- Performance benchmarks"

# Then implement to pass tests:
"Now implement .detectCompositionalBias() to pass all these tests"
```

### B. Example Data Generation
```markdown
"Create realistic example datasets for DAAadvisor that cover:
1. High sparsity 16S data (90%+ zeros)
2. Moderate sparsity metagenomic data
3. Extreme sparsity viral data
4. Data with strong compositional bias
Include metadata with batch effects"
```

## 7. Documentation with Claude

### A. Vignette Development
```markdown
"Help me write a vignette section on 'Understanding Method Recommendations'.
It should:
- Explain how the selection algorithm works
- Show real examples with different data types
- Include interpretation guidelines
- Be accessible to non-statisticians"
```

### B. Function Documentation
```markdown
"Review and improve this roxygen2 documentation:
[paste current docs]

Ensure it includes:
- Clear description
- All parameters explained
- Return value structure
- Examples that actually run
- Links to related functions"
```

## 8. Performance Optimization

### A. Profiling Requests
```markdown
"This function is slow with large datasets:
[paste code]

Please:
1. Identify bottlenecks
2. Suggest vectorized alternatives
3. Consider Rcpp if needed
4. Maintain readability"
```

### B. Memory Efficiency
```markdown
"The data profiler uses too much memory with 50k features.
How can I:
1. Process in chunks
2. Use sparse matrices where appropriate
3. Clean up intermediate objects
4. Add progress bars for long operations"
```

## 9. Package Finalization

### A. CRAN Preparation
```markdown
"Help me prepare for CRAN submission:
1. Check my DESCRIPTION file
2. Ensure all examples run < 5 seconds
3. Add \dontrun{} where needed
4. Fix any NOTEs from R CMD check"
```

### B. User Interface Polish
```markdown
"Make the package output more user-friendly:
1. Improve print methods with better formatting
2. Add color to console output (with cli package)
3. Create informative warning/error messages
4. Add progress indicators"
```

## 10. Continuous Improvement

### A. Feature Requests Template
```markdown
"Users requested support for longitudinal data.
How can I extend the current framework to:
1. Handle repeated measures
2. Detect temporal patterns
3. Adjust method recommendations
Without breaking existing functionality"
```

### B. Bug Fix Template
```markdown
"Bug report: [description]
Reproducible example: [code]
Expected behavior: [description]
Actual behavior: [description]

Please help me:
1. Diagnose the issue
2. Write a failing test
3. Implement the fix
4. Ensure no regression"
```

## Key Success Factors

1. **Maintain Context**: Start each session with clear context
2. **Incremental Progress**: Build and test small pieces
3. **Code Review**: Have Claude review each major component
4. **Test Everything**: Write tests before/after implementation
5. **Document as You Go**: Don't leave docs for later
6. **Version Control**: Commit working versions frequently
7. **Real Data Testing**: Use actual microbiome datasets early

## Example Session Starter

```markdown
# Session Template
"I'm continuing development of DAAadvisor, an R package for intelligent differential abundance analysis.

Progress so far:
- ✓ Package structure created
- ✓ Data profiler 80% complete
- ⚠ Method selector needs scoring implementation
- ⚡ Currently working on: [specific task]

Package philosophy:
- Intelligent method selection based on data characteristics
- Better than DAtest because we analyze actual data, not just spike-ins
- User-friendly with clear explanations

Today's goal: [specific deliverable]

Here's my current code: [paste relevant section]
Here's the design we're following: [paste relevant design section]

Please help me: [specific request]"
```

This strategy ensures efficient, high-quality implementation while maintaining consistency with the original design vision.