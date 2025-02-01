#pragma once

#include <cstdlib>
#include <memory>
#include <vector>
#include <complex>
#include <cmath>
#include <limits>

namespace bitlush
{
	using f64 = double;
	using c64 = std::complex<f64>;

	inline f64 madsen_l2_squared(const c64& c)
	{
		return c.real() * c.real() + c.imag() * c.imag();
	}

	inline f64 madsen_l2_squared(f64 v)
	{
		return v * v;
	}

	inline f64 madsen_l2(const c64& c)
	{
		return std::sqrt(madsen_l2_squared(c));
	}

	inline f64 madsen_l2(f64 v)
	{
		return std::abs(v);
	}

	inline f64 madsen_l1(const c64& c)
	{
		return std::abs(c.real()) + std::abs(c.imag());
	}

	inline f64 madsen_l1(f64 v)
	{
		return std::abs(v);
	}

	template <typename T>
	class madsen_polynomial final
	{
		T* coefficients_{ nullptr };
		std::size_t degree_{ 0 };

	public:
		madsen_polynomial() = default;

		madsen_polynomial(T* coefficients, std::size_t degree) :
			coefficients_{ coefficients }, degree_{ degree }
		{
		}

		std::size_t degree() const { return degree_; }

		const T& operator[](std::size_t i) const { return coefficients_[i]; }

		void assign(const T* coefficients, std::size_t degree)
		{
			for (std::size_t i = 0; i <= degree; i++)
			{
				coefficients_[i] = coefficients[i];
			}

			degree_ = degree;
		}

		c64 horner(const c64& z) const
		{
			c64 fz = coefficients_[0];

			for (std::size_t i = 1; i <= degree_; i++)
			{
				fz = fz * z + coefficients_[i];
			}

			return fz;
		}

		void deflate(const T& root)
		{
			for (std::size_t k = 1; k < degree_; k++)
			{
				coefficients_[k] += coefficients_[k - 1] * root;
			}

			--degree_;
		}

		void deflate_conjugate_pair(const c64& root)
		{
			T r = -2.0 * root.real();
			T u = madsen_l2_squared(root);

			coefficients_[1] -= r * coefficients_[0];

			for (std::size_t k = 2; k < degree_ - 1; k++)
			{
				coefficients_[k] -= r * coefficients_[k - 1] + u * coefficients_[k - 2];
			}

			degree_ -= 2;
		}

		void differentiate_to(madsen_polynomial& target) const
		{
			for (std::size_t i = 0; i < degree_; i++)
			{
				target.coefficients_[i] = static_cast<f64>(degree_ - i) * coefficients_[i];
			}

			target.degree_ = degree_ - 1;
		}

		void scale()
		{
			f64 max_norm = 0;
			f64 min_norm = std::numeric_limits<f64>::max();

			for (std::size_t i = 0; i <= degree_; i++)
			{
				f64 norm = madsen_l1(coefficients_[i]);

				if (norm > 0)
				{
					max_norm = std::max(norm, max_norm);
					min_norm = std::min(norm, min_norm);
				}
			}

			f64 scale = 1 / std::sqrt(max_norm * min_norm);

			for (std::size_t i = 0; i <= degree_; i++)
			{
				coefficients_[i] *= scale;
			}
		}

		void trim()
		{
			while (degree_ > 0)
			{
				if (madsen_l1(coefficients_[0]) == 0)
				{
					++coefficients_;
					--degree_;
				}
				else if (madsen_l1(coefficients_[degree_]) == 0)
				{
					--degree_;
				}
				else
				{
					break;
				}
			}
		}
	};

	template <typename T>
	class madsen_polynomial_evaluation
	{
		using polynomial = madsen_polynomial<T>;

		const polynomial& p_;
		c64 value_;
		f64 abs_squared_;

	public:
		madsen_polynomial_evaluation(const polynomial& p) :
			p_{ p }
		{
		}

		madsen_polynomial_evaluation(const polynomial& p, const c64& v) :
			p_{ p }, value_{ p.horner(v) }, abs_squared_{ madsen_l2_squared(value_) }
		{
		}

		madsen_polynomial_evaluation(const madsen_polynomial_evaluation&) = default;

		madsen_polynomial_evaluation& operator=(const madsen_polynomial_evaluation& other)
		{
			value_ = other.value_;
			abs_squared_ = other.abs_squared_;

			return *this;
		}

		madsen_polynomial_evaluation& operator=(const c64& v)
		{
			value_ = p_.horner(v);
			abs_squared_ = madsen_l2_squared(value_);

			return *this;
		}

		operator const c64&() const { return value_; }

		const c64& value() const { return value_; }

		const f64& abs_squared() const { return abs_squared_; }
	};

	template <typename T>
	class madsen_root_guess final
	{
		using polynomial = madsen_polynomial<T>;

		const polynomial& p0_;
		const polynomial& p1_;
		c64 dz_;
		c64 z1_;

	public:
		madsen_root_guess(const polynomial& p0, const polynomial& p1) :
			p0_{ p0 }, p1_{ p1 }
		{
		}

		void fast_guess()
		{
			const auto p1z0 = p1_.horner(0);

			if (madsen_l1(p1z0) > 0)
			{
				dz_ = -p0_.horner(0) / p1z0;
			}
			else
			{
				dz_ = 1;
			}

			z1_ = dz_;
		}

		void slow_guess()
		{
			dz_ = initial_dz();
			z1_ = 0.5 * find_min_p_norm() * dz_ / l2_norm(dz_);
		}

		const c64& dz() const { return dz_; }

		const c64& z1() const { return z1_; }

	private:
		c64 initial_dz() const
		{
			const c64 z0 = 0;

			const c64 p1z0 = p1_.horner(z0);

			if (madsen_l1(p1z0) > 0)
			{
				const c64 p0z0 = p0_.horner(z0);

				return -p0z0 / p1z0;
			}
			else
			{
				return 1.0;
			}
		}

		f64 find_min_p_norm() const
		{
			const auto degree = p0_.degree();

			f64 min_ak = std::numeric_limits<f64>::max();

			const f64 a0 = madsen_l2_squared(p0_[degree]);

			for (std::size_t i = 0; i < degree; i++)
			{
				f64 ak = madsen_l2_squared(p0_[i]);

				if (ak != 0)
				{
					ak = std::log(a0 / ak) / static_cast<f64>(2 * (degree - i));

					if (ak < min_ak)
					{
						min_ak = ak;
					}
				}
			}

			return std::exp(min_ak);
		}
	};

	template <typename T>
	class madsen_search_point
	{
		using polynomial = madsen_polynomial<T>;
		using polynomial_evaluation = madsen_polynomial_evaluation<T>;

		const polynomial& p0_;
		c64 value_;
		polynomial_evaluation eval_{ p0_ };

	public:
		madsen_search_point(const polynomial& p0) :
			p0_{ p0 }, value_{ 0 }
		{
		}

		madsen_search_point(const polynomial& p0, const c64& v) :
			p0_{ p0 }, value_{ v }, eval_{ p0, v }
		{
		}

		madsen_search_point(const madsen_search_point&) = default;

		madsen_search_point& operator=(const madsen_search_point& other)
		{
			value_ = other.value_;
			eval_ = other.eval_;

			return *this;
		}

		madsen_search_point& operator=(const c64& v)
		{
			value_ = v;
			eval_ = v;

			return *this;
		}

		madsen_search_point& operator+=(const c64& dz)
		{
			value_ += dz;
			eval_ = value_;

			return *this;
		}

		madsen_search_point operator+(const c64& dz)
		{
			return { p0_, value_ + dz };
		}

		operator const c64&() const { return value_; }

		const c64& value() const { return value_; }

		const polynomial_evaluation& eval() const { return eval_; }
	};

	template <typename T>
	class madsen_root_iterator final
	{
		using search_point = madsen_search_point<T>;
		using polynomial = madsen_polynomial<T>;
		using polynomial_evaluation = madsen_polynomial_evaluation<T>;
		using root_guess = madsen_root_guess<T>;

		static constexpr f64 epsilon = std::numeric_limits<f64>::epsilon();
		static constexpr std::size_t max_iterations = 30;

		c64 theta_ = { 0.8, 0.6 };
		const polynomial& p0_;
		const polynomial& p1_;
		std::size_t iterations_{ 0 };
		search_point z0_{ p0_ };
		search_point z1_{ p0_ };
		search_point zz0_{ p1_ };
		search_point zz1_{ p1_ };
		c64 dz_;
		f64 dz_abs_squared_;
		f64 stopping_criteria_upper_bound_;

	public:
		madsen_root_iterator(const polynomial& p0, const polynomial& p1) :
			p0_{ p0 }, p1_{ p1 }
		{
			root_guess guess{ p0_, p1_ };

			guess.fast_guess();

			z1_ = guess.z1();
			dz_ = guess.dz();
			dz_abs_squared_ = madsen_l2_squared(dz_);

			stopping_criteria_upper_bound_ = stopping_criteria_upper_bound();
		}

		void find_root()
		{
			bool moved = true;
			bool newton = false;

			while (more_iterations_needed())
			{
				if (moved && has_converged()) // convergence check only needed if moved
				{
					break;
				}

				if (newton)
				{
					newton = newton_step();

					if (!newton)
					{
						moved = false; // switch back to search in a new direction
					}
				}
				else
				{
					if (moved)
					{
						tentative_step();
					}
					else
					{
						change_direction();
					}

					moved = search();

					if (moved)
					{
						dz_abs_squared_ = madsen_l2_squared(dz_);

						zz0_ = z0_.value();
						zz1_ = z1_.value();

						newton = use_newton(); // switch to newton
					}
				}
			}
		}

		const c64& root() const { return z1_; }

		const polynomial_evaluation& eval() const { return z1_.eval(); }

	private:
		bool search()
		{
			if (multiplicity_search())
			{
				return true;
			}
			else
			{
				return bisection_search();
			}
		}

		bool bisection_search()
		{
			bool moved = false;

			auto best = z1_ + dz_;

			for (std::size_t p = 1; p <= p0_.degree(); p++)
			{
				const auto halfway = dz_* 0.5;

				const auto bisection = z1_ + halfway;

				if (bisection.eval().abs_squared() < best.eval().abs_squared())
				{
					best = bisection;

					dz_ = halfway;

					moved = true;
				}
				else
				{
					break;
				}
			}

			if (moved)
			{
				z0_ = z1_;
				z1_ = best;
			}

			return moved;
		}

		bool multiplicity_search()
		{
			auto best = z1_ + dz_;

			if (best.eval().abs_squared() < z1_.eval().abs_squared()) //eq. 2.4
			{
				z0_ = z1_;
				z1_ = best;

				for (std::size_t p = 2; p <= p0_.degree(); p++)
				{
					best += dz_;

					if (best.eval().abs_squared() < z1_.eval().abs_squared())
					{
						z1_ = best;
					}
					else
					{
						break;
					}
				}

				dz_ = z1_.value() - z0_.value();

				return true;
			}

			return false;
		}

		void change_direction()
		{
			dz_ *= 0.5 * theta_;
		}

		void tentative_step()
		{
			const c64 f1 = p1_.horner(z1_);

			if (madsen_l1(f1) > 0)
			{
				const c64 f0 = z1_.eval();

				c64 nz = -f0 / f1;

				const f64 r1 = madsen_l2_squared(nz);
				const f64 r0 = madsen_l2_squared(dz_);

				if (r1 > 3.0 * 3.0 * r0)
				{
					nz *= 3.0 * theta_ * std::sqrt(r0 / r1);  // eq. 2.3b
				}

				dz_ = nz;
			}
			else
			{
				change_direction();
			}
		}

		bool newton_step()
		{
			const c64 f1 = zz1_.eval();

			if (madsen_l1(f1) > 0)
			{
				const c64 f0 = z1_.eval();

				const c64 nz = -f0 / f1;

				const auto best = z1_ + nz;

				if (best.eval().abs_squared() < z1_.eval().abs_squared())
				{
					z0_ = z1_;
					z1_ = best;
					dz_ = nz;
					dz_abs_squared_ = madsen_l2_squared(nz);

					zz0_ = zz1_;
					zz1_ = z1_.value();

					return use_newton();
				}
				else
				{
					return false;
				}
			}
			else
			{
				return false;
			}
		}

		bool more_iterations_needed()
		{
			if (iterations_ < max_iterations)
			{
				++iterations_;

				return true;
			}

			return false;
		}

		f64 stopping_criteria_upper_bound() const
		{
			const auto degree = p0_.degree();

			const f64 factor = 16.0 * epsilon * static_cast<f64>(degree);

			return factor * factor * madsen_l2_squared(p0_[degree]);
		}

		bool has_converged() const
		{
			if (madsen_l2_squared(dz_) <= epsilon * epsilon * madsen_l2_squared(z1_))
			{
				return true; // eq. 2.10
			}

			if (z1_.eval().abs_squared() <= stopping_criteria_upper_bound_)
			{
				return true; // eq. 2.11
			}

			return false;
		}

		bool use_newton() const
		{
			const c64 p1z0 = zz0_.eval();
			const c64 p1z1 = zz1_.eval();

			const f64 p0z1_as = z1_.eval().abs_squared();
			const f64 p1z1_as = zz1_.eval().abs_squared();

			return p0z1_as * madsen_l2_squared(p1z0 - p1z1) <= p1z1_as * p1z1_as * madsen_l2_squared(dz_) * 0.25;
		}
	};

	template <typename T>
	class madsen_root_finder;

	template <>
	class madsen_root_finder<f64> final
	{
		using T = f64;

		using root_iterator = madsen_root_iterator<T>;
		using polynomial = madsen_polynomial<T>;

		static constexpr f64 epsilon = std::numeric_limits<f64>::epsilon();

		std::vector<T> p0_data_;
		std::vector<T> p1_data_;
		polynomial p0_;
		polynomial p1_;
		std::vector<c64> roots_;

	public:
		const c64& operator[](std::size_t i) const { return roots_[i]; }

		void find_roots(const T* coefficients, std::size_t degree)
		{
			assign(coefficients, degree);

			while (p0_.degree() > 2)
			{
				p0_.scale();
				p0_.differentiate_to(p1_);

				root_iterator iterator{ p0_, p1_ };

				iterator.find_root();

				const auto root = iterator.root();

				const bool real = is_real_root(root, std::sqrt(iterator.eval().abs_squared()));

				if (real)
				{
					roots_[p0_.degree() - 1] = root.real();

					p0_.deflate(root.real());
				}
				else
				{
					roots_[p0_.degree() - 1] = root;
					roots_[p0_.degree() - 2] = std::conj(root);

					p0_.deflate_conjugate_pair(root);
				}
			}

			if (p0_.degree() == 2)
			{
				solve_quadratic();
			}
			else if (p0_.degree() == 1)
			{
				solve_linear();
			}
		}

	private:
		void assign(const T* coefficients, std::size_t degree)
		{
			p0_data_.resize(degree + 1);
			p1_data_.resize(degree);

			p0_ = { p0_data_.data(), p0_data_.size() };
			p1_ = { p1_data_.data(), p1_data_.size() };

			p0_.assign(coefficients, degree);
			p0_.trim();

			roots_.resize(degree);

			for (std::size_t i = 0; i < p0_.degree(); i++)
			{
				roots_[i] = 0;
			}
		}

		bool is_real_root(const c64& root, f64 abs_f) const
		{
			const f64 x = root.real();
			const f64 abs_x = std::abs(x);

			f64 s = p0_[0];
			f64 g = 0;
			f64 abs_s = std::abs(s);

			for (std::size_t i = 1; i <= p0_.degree(); i++)
			{
				const f64 ss = x * s + p0_[i];
				const f64 abs_ss = std::abs(ss);
				const f64 gg = abs_x * (g + abs_s) + abs_ss;

				g = gg;
				s = ss;
				abs_s = abs_ss;
			}

			return abs_s < 2 * epsilon * g + abs_f;
		}

		void solve_linear() // a*x + b = 0
		{
			roots_[0] = -p0_[1] / p0_[0];
		}

		void solve_quadratic() // a*x^2 + b*x + c = 0
		{
			const f64 a = p0_[0];
			const f64 b = p0_[1];
			const f64 c = p0_[2];

			const c64 discriminant = std::sqrt(c64(b * b - 4 * a * c));

			const c64 r = b * discriminant;

			if (r.real() >= 0)
			{
				const c64 q = -0.5 * (b + discriminant);

				roots_[0] = q / a;
				roots_[1] = c / q;
			}
			else
			{
				const c64 q = -0.5 * (b - discriminant);

				roots_[0] = q / a;
				roots_[1] = c / q;
			}
		}
	};

	template <>
	class madsen_root_finder<c64> final
	{
		using T = c64;

		using root_iterator = madsen_root_iterator<T>;
		using polynomial = madsen_polynomial<T>;

		std::vector<T> p0_data_;
		std::vector<T> p1_data_;
		polynomial p0_;
		polynomial p1_;
		std::vector<c64> roots_;

	public:
		const c64& operator[](std::size_t i) const { return roots_[i]; }

		void find_roots(const T* coefficients, std::size_t degree)
		{
			assign(coefficients, degree);

			while (p0_.degree() > 2)
			{
				p0_.scale();
				p0_.differentiate_to(p1_);

				root_iterator iterator{ p0_, p1_ };

				iterator.find_root();

				const auto root = iterator.root();

				roots_[p0_.degree() - 1] = root;

				p0_.deflate(root);
			}

			if (p0_.degree() == 2)
			{
				solve_quadratic();
			}
			else if (p0_.degree() == 1)
			{
				solve_linear();
			}
		}

	private:
		void assign(const T* coefficients, std::size_t degree)
		{
			p0_data_.resize(degree + 1);
			p1_data_.resize(degree);

			p0_ = { p0_data_.data(), p0_data_.size() };
			p1_ = { p1_data_.data(), p1_data_.size() };

			p0_.assign(coefficients, degree);
			p0_.trim();

			roots_.resize(degree);

			for (std::size_t i = 0; i < p0_.degree(); i++)
			{
				roots_[i] = 0;
			}
		}

		void solve_linear() // a*x + b = 0
		{
			roots_[0] = -p0_[1] / p0_[0];
		}

		void solve_quadratic() // a*x^2 + b*x + c = 0
		{
			const c64 a = p0_[0];
			const c64 b = p0_[1];
			const c64 c = p0_[2];

			const c64 discriminant = std::sqrt(b * b - 4.0 * a * c);

			const c64 r = std::conj(b) * discriminant;

			if (r.real() > 0)
			{
				const c64 q = -0.5 * (b + discriminant);

				roots_[0] = q / a;
				roots_[1] = c / q;
			}
			else
			{
				const c64 q = -0.5 * (b - discriminant);

				roots_[0] = q / a;
				roots_[1] = c / q;
			}
		}
	};
}

