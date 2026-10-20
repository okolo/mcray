#ifndef MCRAY_CRBEAM_SHA256_H
#define MCRAY_CRBEAM_SHA256_H

#include <string>

// SHA-256 file fingerprint for reproducibility and manifest verification.
// Uses OpenSSL Crypto EVP; throws on read or digest failure.
std::string Sha256File(const std::string &path);

#endif
