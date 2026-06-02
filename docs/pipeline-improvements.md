# Pipeline Improvements Documentation

This document describes five key improvements made to the nf-codonyat pipeline for enhanced reliability, usability, and maintainability.

## Features

### 1. Retry Strategy

The pipeline now implements an intelligent retry mechanism for robust handling of transient failures:

- **Automatic retries** with exponential backoff
- **Configurable retry count** (default: 2 attempts)
- **Task-level granularity** for precise failure handling
- **Error classification** to distinguish transient from permanent failures

**Benefits:**
- Reduces false negatives from temporary system issues
- Improves pipeline success rate without code changes
- Transparent retry behavior with clear logging

### 2. Input Validation

Early validation prevents invalid runs and provides clear feedback:

- **Comprehensive input checking** before pipeline execution
- **Format validation** for sequence data and metadata
- **Range checking** for numerical parameters
- **Dependency validation** between related parameters

**Benefits:**
- Catches configuration errors early
- Reduces debugging time for users
- Prevents wasted compute resources on invalid inputs

### 3. Error Handling

Improved error messages and recovery mechanisms:

- **Structured error messages** with context and suggestions
- **Graceful failure modes** that preserve intermediate results
- **Error categorization** (input validation, resource, system, external)
- **Clear recovery paths** for common failure scenarios

**Benefits:**
- Users understand what went wrong and how to fix it
- Pipeline operators can diagnose issues faster
- Enables better monitoring and alerting

### 4. Conda Profile

Native conda support for dependency management:

- **Conda-based environment** alongside Docker
- **Package resolution** via conda channels
- **Environment reproducibility** with conda-lock files
- **Local development** without Docker overhead

**Benefits:**
- Flexible deployment options
- Better compatibility with HPC environments
- Simplified local testing and development

### 5. Local Package Support

Use unpublished packages from local directories:

- **Local directory references** in configuration
- **Development workflow** integration
- **Path-based package resolution** for testing
- **Seamless switching** between local and published versions

**Benefits:**
- Accelerates development cycles
- Enables integration testing of components
- Simplifies collaborative development

## Usage Examples

### Docker Profile

Standard production deployment using Docker:

```bash
# Basic execution
nextflow run main.nf -profile docker

# With custom resources
nextflow run main.nf -profile docker \
  --max_cpus 8 \
  --max_memory 32.GB

# With retry configuration
nextflow run main.nf -profile docker \
  --retry_attempt_count 3
```

### Conda Profile

Local development or HPC deployment:

```bash
# Using conda environment
nextflow run main.nf -profile conda

# With conda-specific options
nextflow run main.nf -profile conda \
  --conda_cache_dir /shared/conda/cache

# Combined with local packages
nextflow run main.nf -profile conda \
  --local_packages_dir ./local_packages
```

### Local Package Development

Testing unpublished components:

```bash
# Reference local package
nextflow run main.nf -profile docker \
  --local_packages_dir ./aa_caller_dev

# Multiple local packages
nextflow run main.nf -profile conda \
  --local_packages_dir "/path/to/pkg1:/path/to/pkg2"
```

## Common Error Messages and Solutions

| Error Message | Cause | Solution |
|---------------|-------|----------|
| `Input validation failed: missing required field 'samples'` | Configuration missing required parameter | Check `nextflow_schema.json` for required fields and provide value |
| `Sequence format invalid: expected FASTA` | Input file format mismatch | Verify input file format matches specified type |
| `Resource request exceeds limits: requested 64GB, max 32GB` | Memory request too high | Increase `--max_memory` or reduce per-task memory |
| `Task failed after 3 retry attempts` | Persistent task failure | Check logs in `work/` directory and investigate task-specific error |
| `Package not found in conda channels` | Missing dependency | Update conda channels or check spelling |
| `Local package directory not found: /path/to/pkg` | Invalid path reference | Verify path exists and is accessible |
| `Docker image pull failed` | Network or registry issue | Check Docker daemon and registry connectivity |
| `No space left on device` | Insufficient disk space | Clear temporary files or increase available space |

## Testing

### Validation Testing

Test input validation with various scenarios:

```bash
# Test with minimal valid config
nextflow run main.nf -profile docker -c tests/config/minimal.config

# Test with invalid inputs
nextflow run main.nf -profile docker -c tests/config/invalid.config

# Verify error messages
nextflow run main.nf -profile docker -c tests/config/missing_param.config
```

### Retry Testing

Verify retry mechanism:

```bash
# Simulate transient failure
nextflow run main.nf -profile docker \
  -c tests/config/retry_test.config \
  --retry_attempt_count 3

# Check retry logs in work/
find work -name "*.log" | xargs grep "retry"
```

### Conda Testing

Test conda profile:

```bash
# Create conda environment
nextflow run main.nf -profile conda --help

# Verify environment creation
ls -la ~/.nextflow/conda/

# Test package resolution
nextflow run main.nf -profile conda -c tests/config/conda_test.config
```

### Local Package Testing

Verify local package resolution:

```bash
# Copy package to local directory
cp -r aa_caller local_packages/

# Run with local package
nextflow run main.nf -profile conda \
  --local_packages_dir ./local_packages \
  -c tests/config/local_pkg_test.config

# Verify package used
grep "local_packages" work/*/command
```

## Files Modified

| File | Changes |
|------|---------|
| `main.nf` | Input validation block, error handling, retry logic |
| `nextflow.config` | Retry parameters, conda profile, local package paths |
| `conf/base.config` | Task retry configuration, error strategy |
| `conf/profiles.config` | Docker and conda profile definitions |
| `lib/validators.nf` | Input validation functions |
| `lib/error_handlers.nf` | Error handling and recovery logic |

## Design Document Reference

For detailed architectural decisions and design rationale, see:
- [Design Document](./design-document.md) - Complete system design and implementation details
- [Retry Strategy Design](./design-document.md#retry-strategy) - Retry mechanism architecture
- [Validation Framework](./design-document.md#input-validation) - Input validation approach
- [Error Handling Architecture](./design-document.md#error-handling) - Error categorization and responses
- [Profile Configuration](./design-document.md#profiles) - Profile definitions and usage

## Additional Resources

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow Conda Integration](https://www.nextflow.io/docs/latest/conda.html)
- [Error Handling Best Practices](https://www.nextflow.io/docs/latest/process.html#errorstrategy)
