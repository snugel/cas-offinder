#pragma once
//range iterator for use in for loop
//an "each" iterator that allows iteration with iterators
//an "each" iterator that allows iteration with two iterators
//an "over" iterator that allow iteration with two containers
/*
Example:
for(int x:range(100){
	//do something 100 times, x going from 0-99
}
for(int x:range(10,50)){
	//do something with x going from 10-49
}
for(int x:range(50,10,-1){
	//do something iterating backwards from 50 to 11
}
for(int x:range(50,100,5){
	//do something iterating backwards from 50 to 95 in increments of 5
}
for(pair<double *,int *> x : over(vector<double>(10),array<int,20>())){
	//iterates x 10 times through the first 10 items of the vector and the first 10 items of the array
}
vector<double> Vec(100,5);//fill Vec of size 100 with the value 5
for(double & x : each(Vec.rbegin(),Vec.rend())){
	//iterates x 100 times through Vec backwards
}
for(pair<double *,double *> x : each(Vec.begin(),Vec.end(),Vec.rbegin(),Vec.rend())){
	//iterates x.first 100 times through Vec forwards and x.second 100 times though Vec backwards
}
*/
#include <utility>
namespace RangeStuff{
	template<typename Type>
	struct RangeIncr{
		struct RangeIteratorForward{
			Type x;
			bool operator != (RangeIteratorForward & Other){
				return x < Other.x;
			}
			void operator ++ (){
				x++;
			}
			Type & operator *(){
				return x;
			}
		};
		RangeIteratorForward Start;
		RangeIteratorForward End;
		RangeIncr(Type InStart, Type InEnd){
			Start = RangeIteratorForward{ InStart };
			End = RangeIteratorForward{InEnd};
		}
		RangeIteratorForward & begin(){
			return Start;
		}
		RangeIteratorForward & end(){
			return End;
		}
	};
	template<typename Type>
	struct RangeIterate{
		struct RangeIterator{
			Type x;
			Type incr;
			bool operator != (RangeIterator & Other){
				return (incr < 0) ? (Other.x < x) : (x < Other.x);//incr == 0 is an exception that is handled
			}
			void operator ++ (){
				(x += incr);
			}
			Type & operator * (){
				return x;
			}
		};
		RangeIterator Start, End;
		RangeIterate(Type InStart, Type InEnd, Type InIncr){
			Start = RangeIterator{ InStart, InIncr };
			End = RangeIterator{ InEnd, InIncr };
		}
		RangeIterator & begin(){
			return Start;
		}
		RangeIterator & end(){
			return End;
		}
	};
	template<typename ItTy>
	class EachIterate{
	public:
		ItTy Begin;
		ItTy End;
		EachIterate(ItTy InBegin, ItTy InEnd){
			Begin = InBegin;
			End = InEnd;
		}
		ItTy & begin(){
			return Begin;
		}
		ItTy & end(){
			return End;
		}
	};
	template<typename ItTy1,typename ItTy2>
	class DuoEachIterate{
	public:
		class DuoIterate{
		public:
			ItTy1 It1;
			ItTy2 It2;
			DuoIterate(ItTy1 & InIt1, ItTy2 & InIt2){
				It1 = InIt1;
				It2 = InIt2;
			}
			DuoIterate() = default;
			bool operator !=(DuoIterate & Other){
				return It1 != Other.It1 && It2 != Other.It2;
			}
			auto operator * ()->decltype(std::make_pair(&(*It1), &(*It2))){
				return std::make_pair(&(*It1),&(*It2));
			}
			void operator ++ (){
				++It1;
				++It2;
			}
		};
		DuoIterate Begin;
		DuoIterate End;
		DuoEachIterate(ItTy1 InBegin1, ItTy1 InEnd1, ItTy2 InBegin2, ItTy2 InEnd2){
			Begin = DuoIterate(InBegin1, InBegin2);
			End = DuoIterate(InEnd1, InEnd2);
		}
		DuoIterate & begin(){
			return Begin;
		}
		DuoIterate & end(){
			return End;
		}
	};
}
template<typename StartType,typename EndType,typename IncrType>
inline auto range(StartType start,EndType end,IncrType incr)->RangeStuff::RangeIterate<long long>{
	return RangeStuff::RangeIterate<long long>(start,end,incr);
}
template<typename StartType,typename EndType>
inline auto range(StartType start,EndType end) ->RangeStuff::RangeIncr<long long>{
	return RangeStuff::RangeIncr<long long>(start,end);
}
template<typename EndType>
inline RangeStuff::RangeIncr<long long> range(EndType end){
	return RangeStuff::RangeIncr<long long>(0,end);
}
template<typename ItTy>
inline RangeStuff::EachIterate<ItTy> each(ItTy InBegin, ItTy InEnd){
	return RangeStuff::EachIterate<ItTy>(InBegin, InEnd);
}
template<typename ItTy1,typename ItTy2>
inline RangeStuff::DuoEachIterate<ItTy1,ItTy2> each(ItTy1 InBegin1, ItTy1 InEnd1,ItTy2 InBegin2, ItTy2 InEnd2){
	return RangeStuff::DuoEachIterate<ItTy1,ItTy2>(InBegin1, InEnd1,InBegin2, InEnd2);
}
template<typename Cont1,typename Cont2>
inline auto over(Cont1 Container1, Cont2 Container2)->RangeStuff::DuoEachIterate<decltype(Container1.begin()),decltype(Container2.begin())>{
	return RangeStuff::DuoEachIterate<decltype(Container1.begin()),decltype(Container2.begin())>(Container1.begin(), Container1.end(),Container2.begin(), Container2.end());
}
