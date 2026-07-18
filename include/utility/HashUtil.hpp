//
// Created by CanhDo on 2024/12/17.
//

#ifndef HASHUTIL_HPP
#define HASHUTIL_HPP

class HashUtil {
public:
    static std::size_t combinedHash(std::size_t hash1, std::size_t hash2) {
        return hash1 ^ (hash2 + 0x9e3779b97f4a7c15ULL + (hash1 << 6) + (hash1 >> 2));
    }
};

#endif //HASHUTIL_HPP
