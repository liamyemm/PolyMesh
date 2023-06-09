#include <thread>
#include <functional>

#ifndef _PARALLEL_FOR_HPP
#define _PARALLEL_FOR_HPP

/// Parallel for
void parallel_for(const unsigned nb_elements, const std::function<void(size_t start, size_t end)> &functor, bool use_threads = true);

#endif