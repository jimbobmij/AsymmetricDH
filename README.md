# Asymmetric Diffie-Hellman Key Agreement Protocol
This repository contains Java implementations of an **asymmetric variant** of the Diffie-Hellman key agreement protocol, designed to reduce the computational burden on one party (typically Alice) by leveraging parallel computing and short exponents. The implementation supports performance measurements and is suitable for experimental analysis in asymmetric computational environments.

## Overview

The protocol is designed for scenarios where the computational capabilities of the communicating parties are significantly unbalanced. This asymmetric structure enables efficient key agreement by offloading most of the computational work to the more powerful device (e.g., Bob), while maintaining the security properties of traditional D-H protocols.

## Files

### `Asymmetric_DH_together.java`

This implementation simulates both **Alice and Bob** executing the asymmetric D-H key agreement on a single machine. It measures:
- Total computation time for key generation and shared secret key (SSK) computation.
- CPU load and number of cores used.
- Average thread execution time per parallel exponentiation.

This version is suitable for evaluating total performance including Bob's public key generation and SSK computation.

### `Asymmetric_DH_together_Alice_only.java`

This implementation focuses on **Alice’s side only**, assuming Bob's public keys are precomputed and available. It:
- Measures Alice's public key and SSK generation times under various configurations.
- Allows input of different values of `d` to simulate varying degrees of parallelism.
- Estimates CPU utilization and per-thread execution time.

This file is useful for analyzing how parallelism and exponent reduction affect Alice's performance in practice.
