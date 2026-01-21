---
name: code-test
description: Write unit tests for code to test it.
---

# Code Test Skill

When testing code, follow these steps:

## Test checklist

* The usage of a prelude file is prohibited.

* All tests for public methods are in the crate_name/tests folder. 

* The tests folder replicates the exact src folder structure, for for example:

tests/errors/mod.r  # contains tests for each error type in a separate file
tests/traits/mod.rs # Optional contains tests for each  trait in a separate file
tests/type/mod.rs # contains tests for each  type in a separate file

* Every single test files must be registered to the correspoding mod file and that module must be registered with its higher up module.

* Ensure the corrext #[cfg(test)] annotation is set for each registeres test file.

* Test files replicate the source file name with an appended _tests. For example,
a source file

`src/errors/normal_error/normal_error.rs`

is matched with the test file under the tests folder:

`test/errors/normal_error/normal_error_tests.rs`

Shared utils used for testing are actually stored in the src tree under:

`src/utils_tests/mod.rs # contains utils`

* Aim for 100% test coverage of all source code files, all lines of code, all code branches and all error cases.
