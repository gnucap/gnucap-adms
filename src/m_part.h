
/* This file is part of "Gnucap-adms"
 *
 * (c) 2015 Felix Salfelder
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * partitions
 */

template<unsigned U>
class PARTITION{
	public:
		PARTITION()
		{
			for(unsigned i=0; i<U; ++i){
				_p[i].first = _p[i].second = i;
			}
		}
		unsigned identify(unsigned a, unsigned b)
		{
			assert(a<U);
			assert(b<U);
			if(cycle(a)==cycle(b)){ untested();
				return cycle(a);
			}else if(cycle(a) < cycle(b)){
				set_cycle(b, cycle(a));
			}else{
				set_cycle(a, cycle(b));
			}
			std::swap(_p[a].first, _p[b].first);
			return cycle(a);
		}
		unsigned cycle(unsigned a) const{
			return _p[a].second;
		}
	private:
		unsigned& _cycle(unsigned a){
			return _p[a].second;
		}
		void set_cycle(unsigned a, unsigned to)
		{ untested();
			assert(a<U);
			assert(to<U);
			unsigned s=a;
			_cycle(s) = to;
			while(_p[s].first != a) { untested();
				s = _p[s].first;
				_cycle(s) = to;
			}
		}
		std::pair<unsigned,unsigned> _p[U];
};
