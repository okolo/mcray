#include "Sha256.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <openssl/evp.h>
#include <sstream>
#include <stdexcept>

std::string Sha256File(const std::string &path) {
  std::ifstream file(path, std::ios::binary);
  if (!file)
    throw std::runtime_error("Cannot hash " + path);
  // EVP is supported by both OpenSSL 1.1.1 and 3.x; do not use deprecated
  // low-level SHA256_Init/Update/Final entry points.
  const auto context = std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)>(
      EVP_MD_CTX_new(), EVP_MD_CTX_free);
  if (!context || EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) != 1)
    throw std::runtime_error("SHA-256 initialization failed");
  std::array<char, 65536> buffer;
  while (file.read(buffer.data(), buffer.size()) || file.gcount())
    if (EVP_DigestUpdate(context.get(), buffer.data(),
                         static_cast<std::size_t>(file.gcount())) != 1)
      throw std::runtime_error("SHA-256 update failed");
  if (!file.eof() || file.bad())
    throw std::runtime_error("Cannot read " + path);
  std::array<unsigned char, EVP_MAX_MD_SIZE> digest;
  unsigned length = 0;
  if (EVP_DigestFinal_ex(context.get(), digest.data(), &length) != 1 ||
      length != 32)
    throw std::runtime_error("SHA-256 finalization failed");
  std::ostringstream out;
  for (unsigned i = 0; i < length; ++i)
    out << std::hex << std::setw(2) << std::setfill('0')
        << static_cast<unsigned>(digest[i]);
  return out.str();
}

std::string Sha256Bytes(const std::string& bytes) {
  std::array<unsigned char, EVP_MAX_MD_SIZE> digest;
  unsigned length = 0;
  if (EVP_Digest(bytes.data(), bytes.size(), digest.data(), &length, EVP_sha256(), nullptr) != 1 || length != 32)
    throw std::runtime_error("SHA-256 byte digest failed");
  std::ostringstream out;
  for (unsigned i=0; i<length; ++i)
    out << std::hex << std::setw(2) << std::setfill('0') << static_cast<unsigned>(digest[i]);
  return out.str();
}
