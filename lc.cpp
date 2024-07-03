// lc.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// date: 2020-02-20
// last modification: 2020-04-16
//

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "ObjBase.h"
#include <vector>
#include "limits.h"
#include <random>
#include <cmath>
#include <queue>

/* 2020-02-20: 开始写label cut程序了！*/


// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件

using namespace std;


/////////////////////////////////////////////////////////////////////
// 类型定义（1）
/////////////////////////////////////////////////////////////////////

typedef set<int> int_set_t;	// 整数的集合
typedef vector<int> int_vector_t;	// 整数向量
typedef vector<bool> bool_vector_t;	// 布尔向量
typedef vector< pair<int, int> > int_pair_vector_t; // 整数对向量

typedef int vertex_t;	// 顶点
typedef string label_t;	// 标签
typedef int weight_t;	// 权重

typedef set<vertex_t> vertex_set_t;	// 顶点集
typedef vector<vertex_t> vertex_vector_t;	// 顶点向量
typedef set<label_t> label_set_t;	// 标签集
typedef vector<label_t> label_vector_t;	// 标签向量
typedef set<weight_t> weight_set_t;	// 权重的集合


/////////////////////////////////////////////////////////////////////
// 全局变量
/////////////////////////////////////////////////////////////////////

// 使用缺省的随机引擎作为随机源
std::default_random_engine dre;


/////////////////////////////////////////////////////////////////////
// 全局函数（1）
/////////////////////////////////////////////////////////////////////

string GenerateGuid();


/////////////////////////////////////////////////////////////////////
// 类定义（1）
/////////////////////////////////////////////////////////////////////

// s和t
// 各种类型的标签s-t割问题都有s和t
// 路的源和目标也是s和t
// 所以把s和t单独作成一个类
class st_t
{
public:
	vertex_t s;	// 源
	vertex_t t;	// 目标

public:
	st_t() :s(0), t(0) {}
};

// 无向图上的边
class edge_t
{
private:
	// 对于无向图，i, j表示边的两个端点。
	// 注意：edge_t不保证i总是小于等于j
	vertex_t i, j;	

public:
	edge_t():i(0),j(0){}
	edge_t(vertex_t ii, vertex_t jj) : i(ii), j(jj) {}
	// 拷贝构造函数
	edge_t(const edge_t& a) : i(a.i), j(a.j){}
	// 拷贝赋值
	edge_t& operator=(const edge_t&);
	vertex_t get_i() const { return i; }
	vertex_t get_j() const { return j; }
	// 将edge_t对象放入自动排序的容器中时，需要提供比较操作符“<”。
	bool operator< (const edge_t& rhs) const;
	bool operator> (const edge_t& rhs) const;
	bool operator== (const edge_t& rhs) const;
};

// 多重图上的边
// 可用于表示有向边，也可以用于表示无向边
class medge_t :public edge_t{
private:
	string id;	// 多重图上的重边，用私有的id进行区分

public:
	medge_t() : edge_t(0, 0) {}
	medge_t(vertex_t ii, vertex_t jj) : edge_t(ii, jj) { id = GenerateGuid(); }
	// 拷贝构造函数
	medge_t(const medge_t& a) : edge_t(a), id(a.id){}
	// 拷贝赋值
	medge_t& operator=(const medge_t&);
	// 将medge_t对象放入自动排序的容器中时，需要提供比较操作符“<”。
	bool operator< (const medge_t &rhs) const;
};

typedef set<edge_t> edge_set_t;	// 边的集合
typedef vector<edge_t> edge_vector_t;	// 边向量
typedef set<medge_t> medge_set_t;	// 多重边的集合
typedef set<medge_t> medge_vector_t;	// 多重边向量
typedef list<medge_t> medge_list_t;	// 多重边链表

// mpatht是多重图上的路
class mpath_t : public st_t {
public:
	medge_list_t e_list;	// 将路径定义为边的序列
};


/////////////////////////////////////////////////////////////////////
// 类型定义（2）
/////////////////////////////////////////////////////////////////////

//typedef edge_vector_t path_t;	// 将路径定义为边的序列
//typedef medge_vector_t mpath_t;	// 将路径定义为边的序列

// 定义多重图的邻接表
// 邻接表的结构为：<顶点 - 和此顶点关联的所有边的集合>，...
// 即，每个顶点和一个边集构成一对。边集中存放的是以该顶点为端点的所有边。
// 所有这样的<顶点, 边集>对构成一个map，即邻接表。
typedef map<vertex_t, edge_set_t> vertex_edge_set_map_t;

// 定义多重图的邻接表
// 邻接表的结构为：<顶点 - 和此顶点关联的所有边的集合>，...
// 即，每个顶点和一个边集构成一对。边集中存放的是以该顶点为端点的所有边。
// 所有这样的<顶点, 边集>对构成一个map，即邻接表。
typedef map<vertex_t, medge_set_t> vertex_medge_set_map_t;	

// 从多重边到标签的映射
typedef map<edge_t, label_t> edge_label_map_t;
// 从多重边到标签的映射
typedef map<medge_t, label_t> medge_label_map_t;	
// 从多重边到权重的映射
typedef map<medge_t, weight_t> medge_weight_map_t;
// 从标签到权重的映射
typedef map<label_t, weight_t> label_weight_map_t;


/////////////////////////////////////////////////////////////////////
// 类定义（2）
/////////////////////////////////////////////////////////////////////

// 无向图
class graph_t
{
public:
	edge_set_t E;
	vertex_edge_set_map_t ve_map;
	
public:
	int vertex_number() const;
	int edge_number() const;
};

// 无向多重图
class multigraph_t 
{
public:
	// 重边的集合E
	// 注意：E仍是一个set，而不是一个multiset
	// 在实现重边（medge_t）时，重边（有id）之间仍是可区分的
	medge_set_t E;
	// ve_map：邻接表，记录每个顶点关联的边
	// ve_map是冗余结构，每条边e在其两个端点处各出现一次
	vertex_medge_set_map_t ve_map;

public:
	int vertex_number() const;
	int edge_number() const;
	int edge_multiplicity() const;
	int edge_number_between_vertices(vertex_t u, vertex_t v) const;
	medge_set_t medges_between_vertices(vertex_t u, vertex_t v) const;
};

// 边上带权重的无向多重图
class weighted_multigraph_t: public multigraph_t 
{
public:
	medge_weight_map_t w;	// 边集到权重的映射

	// 定义一个以父类对象为参数的拷贝构造函数，
	// 以允许将父类（multigraph_t）对象赋值给子类（weighted_multigraph_t）对象
	weighted_multigraph_t(const multigraph_t& rhs)
	{
		this->E = rhs.E;
		this->ve_map = rhs.ve_map;
	}
};

// 标签s-t割问题，标签上带权重
class weighted_labelstcut_t : public st_t
{
public:
	graph_t G;	// 多重图
	label_set_t L;	// 标签集合

	edge_label_map_t ell;	// 边集到标签集的映射
	label_weight_map_t w;	// 标签集到权重的映射

public:
	weight_t max_weight();
};

// 多重图上的标签s-t割问题，标签上没有权重，带禁止标签
class labelstcut_m_f_t : public st_t
{
public:
	multigraph_t G;	// 多重图
	label_set_t A;	// 在解中允许使用的标签的集合
	label_set_t B;	// 在解中不允许使用的标签的集合

	medge_label_map_t ell;	// 边集到标签集的映射

public:
	label_set_t labels(medge_set_t E);
};

// 多重图上的标签s-t割问题，标签上有权重，带禁止标签
class labelstcut_w_m_f_t : public labelstcut_m_f_t
{
public:
	label_weight_map_t w;	// 标签集到权重的映射
};

// 枚举器
class enumerator_t
{
private:
	int n;	// 元素的数目
	int_set_t current_enum;

public:
	enumerator_t() {
		n = 0;
		current_enum.clear();
	}
	enumerator_t(int nn) : n(nn) { current_enum.clear(); }
	void initialize(int n);
	int_set_t next_enumeration();
};


/////////////////////////////////////////////////////////////////////
// 全局函数（2）
/////////////////////////////////////////////////////////////////////

// 得到多重图的顶点数
int vertex_number(const multigraph_t& G);

// 得到多重图的边数
int edge_number(const multigraph_t& G);

// 得到多重图的边的最大重数
int edge_multiplicity(const multigraph_t& G);

// 将标签集转换为标签向量
label_vector_t label_set_to_label_vector(const label_set_t& L);

// 按G(n, p)模型产生一个随机图
graph_t generate_a_random_graph_Gnp
(
	int n, // n -- 顶点数目
	double p // p -- 两点之间有边的概率。
);

// 按配置模型产生一个随机图
graph_t generate_a_random_graph_CM
(
	int n // n -- 顶点数目
);

// 产生一个随机的多重图
multigraph_t generate_a_random_multigraph
(
	int n,	// 顶点数目
	double p, // 两点之间有边的概率
	int mu	// 两点之间产生边时，重复的次数
);

// 产生一个随机的边上带权重的多重图
weighted_multigraph_t generate_a_random_weighted_multigraph
(
	int n,	// 顶点数目
	double p, // 两点之间有边的概率
	int mu,	// 两点之间产生边时，重复的次数
	int w_max	// 最大权重
);

// 求G上从s到t的一条最短路
mpath_t shortest_path_w_m
(
	const weighted_multigraph_t& G, 
	vertex_t s, 
	vertex_t t
);

// 求带权多重图G上路p的长度
weight_t path_length
(
	const weighted_multigraph_t& G, 
	const mpath_t& p
);

// 求带权多重图G上从s到t的距离
weight_t dist
(
	const weighted_multigraph_t& G, 
	vertex_t s, 
	vertex_t t
);

// 产生一个随机的weighted_labelstcut_t实例
// 按照ER模型中的G(n, p)模型产生随机图
weighted_labelstcut_t generate_a_random_weighted_labelstcut_instance_Gnp
(
	int n,	// 顶点数目
	int q,	// 标签数目
	double p,	// 两点之间有边的概率
	int w_max	// 标签的最大权重
);

// 产生一个随机的weighted_labelstcut_t实例
// 按照ER模型中的CM模型（配置模型）产生随机图
weighted_labelstcut_t generate_a_random_weighted_labelstcut_instance_CM
(
	int n,	// 顶点数目
	int q,	// 标签数目
	int w_max	// 标签的最大权重
);

// 产生一个随机的labelstcut_w_m_f实例
labelstcut_w_m_f_t generate_a_random_labelstcut_w_m_f_instance
(
	int n,	// 顶点数目
	int q,	// 标签数目
	double p,	// 两点之间有边的概率
	int mu,	// 两个顶点之间产生边时，重复的次数
	int w_max	// 最大权重
);

// 产生q个标签
label_set_t generate_a_label_set(int q);

// 求解（近似算法）多重图上带禁止标签的标签s-t割问题
label_set_t solve_labelstcut_on_multigraphs_with_forbidden_labels
(
	const labelstcut_m_f_t& I
);

// 求解（近似算法）标签上带权重的标签s-t割问题
label_set_t solve_weighted_labelstcut
(
	const weighted_labelstcut_t& I
);

// 将边集E中具有L中标签的边删除
// 注意E是出口参数
void remove_edges_with_specified_labels
(
	edge_set_t& E,
	const edge_label_map_t& ell,
	const label_set_t& L_delete
);


// 将图G中具有L_delete中标签的边删除
// 注意G是出口参数
void remove_edges_with_specified_labels
(
	graph_t& G,
	const edge_label_map_t& ell,
	const label_set_t& L_delete
);

// 将多重图G中具有L_delete中标签的边删除
// 注意G是出口参数
void remove_medges_with_specified_labels
(
	multigraph_t& G,
	const medge_label_map_t& ell,
	const label_set_t& L_delete
);

// 将重边的集合E中具有L_delete中标签的边删除
// 注意E是出口参数
void remove_medges_with_specified_labels
(
	medge_set_t& E,
	const medge_label_map_t& ell,
	const label_set_t& L_delete
);

// 将实例I中的多重图I.G转换为边上带权重的多重图G
// 对I.G中的每一条边，若其标签是允许标签，则该边的权重为1
// 否则（其标签为禁止标签）该边的权重为0
weighted_multigraph_t label_weighted_multigraph_to_edge_weighted_multigraph
(
	const labelstcut_m_f_t& I
);

// 将边集E映射到标签集
label_set_t labels
(
	const medge_label_map_t& ell, 
	const medge_set_t& E
);

// 多重图上带禁止标签的标签s-t割问题的两阶段解法
// 第一阶段，连续去掉最短路
label_set_t stage_1_m_f
(
	weighted_multigraph_t& G,
	vertex_t s,
	vertex_t t,
	const medge_label_map_t& ell,
	const label_set_t& A,	// A是允许标签的集合
	unsigned int OPT
);

// 多重图上带禁止标签的标签s-t割问题的两阶段解法
// 第二阶段，对顶点分层，找一个s-t割
label_set_t stage_2_m_f
(
	weighted_multigraph_t& R,
	vertex_t s,
	vertex_t t,
	const medge_label_map_t& ell
);

// 多重图上带禁止标签的标签s-t割问题的两阶段解法，核心程序
label_set_t core_m_f
(
	weighted_multigraph_t G,
	vertex_t s,
	vertex_t t,
	const medge_label_map_t& ell,
	const label_set_t& A,	// A是允许标签的集合
	unsigned int OPT
);

// 根据I、q_star和W，构造带禁止标签的多重图上的标签s-t割实例
labelstcut_m_f_t labelstcut_weighted_to_mf
(
	const weighted_labelstcut_t& I,
	int q_star,
	weight_t W
);

// 根据ell构造group_ell
label_set_t ell_to_group_ell
(
	label_t ell, 
	weight_t w_ell, 
	int q_star, 
	weight_t W
);

// 识别出L中的允许标签
label_set_t admissible_labels
(
	const label_set_t& L,
	const label_weight_map_t& w, 
	weight_t W
);

// 将group labels的集合转换为labels的集合
label_set_t group_labels_to_labels
(
	const label_set_t& A_tilde_prime,
	const label_set_t& L,
	const label_weight_map_t& w,
	int q_star,
	weight_t W
);

// 使用枚举的方法，求标签带权重的标签st割问题的最优解
label_set_t find_optimal_solution
(
	const weighted_labelstcut_t& I
);

// 判断在G中s和t是否连通
bool connected
(
	const graph_t& G,
	vertex_t s,
	vertex_t t
);

// 从s开始BFS，遇到t则结束
// 返回值：BFS过程中遇到了t则返回真，否则返回假
bool BFS_conditional
(
	const graph_t& G,
	vertex_t s,
	vertex_t t
);

// 计算标签集的权重
weight_t weight
(
	const label_set_t& L, 
	const label_weight_map_t& w
);

// 以G(n, p)模型产生随机实例
// 测试标签带权重的标签s-t割问题的近似解法的性能
void test_Gnp();

// 以配置模型产生随机实例
// 测试标签带权重的标签s-t割问题的近似解法的性能
void test_CM();

// 为L中的每个标签产生一个随机的权重，权重取值为[1, w_max]
label_weight_map_t generate_label_weights
(
	const label_set_t& L,
	weight_t w_max
);

// 为图上的边产生随机的带权重的标签
// 分为如下步骤：
// （1）产生q个标签，并为标签产生随机的权重
// （2）为每条边分配一个随机的标签
void generate_random_weighted_labels_for_edges
(
	weighted_labelstcut_t& I,
	int q,
	weight_t w_max
);

// 产生随机的度序列
int_vector_t generate_a_random_degree_sequence
(
	int n	// 顶点数目
);

// 将d中剩余度大于0的顶点加入到c
// 将剩余度大于0的顶点和它的每个度构成的pair <i, k>加入到u
void translate_remaining_degrees
(
	int n,	// 顶点数目
	const int_vector_t& d,
	int_vector_t& c,	// 出口参数
	int_pair_vector_t& u	// 出口参数
);


/////////////////////////////////////////////////////////////////////
// 各函数的实现代码
/////////////////////////////////////////////////////////////////////


// 测试代码
void test(void)
{
	// 测试代码，产生随机的标签带权重的标签s-t割问题实例
	/*
	weighted_labelstcut_t I = generate_a_random_weighted_labelstcut_instance(
		6,	// 顶点数目
		10,	// 标签数目
		0.6,	// 两点之间有边的概率
		100	// 标签的最大权重
	);
	*/

	// 测试代码，产生随机的边上带权的多重图
	/*
	weighted_multigraph_t G = generate_a_random_weighted_multigraph(
		6,	// 顶点数目
		0.6, // 两点之间有边的概率
		5,	// 两点之间产生边时，重复的次数
		100	// 最大权重
	);
	*/

	// 测试代码，求边上带权的多重图上两点之间的最短路
	/*
	vertex_t s = 0;
	vertex_t t = 2;
	mpath_t p = shortest_path_w_m(G, s, t);
	*/

	// 测试代码，产生随机的标签带权重的标签s-t割问题实例
	/*
	weighted_labelstcut_t I = generate_a_random_weighted_labelstcut_instance(
		6,	// 顶点数目
		10,	// 标签数目
		0.6,	// 两点之间有边的概率
		100	// 标签的最大权重
	);
	*/

	// 测试代码，求解标签带权重的标签s-t割问题
	/*
	weighted_labelstcut_t I = generate_a_random_weighted_labelstcut_instance(
		10,	// 顶点数目
		10,	// 标签数目
		0.4,	// 两点之间有边的概率
		100	// 标签的最大权重
	);

	label_set_t sol = solve_weighted_labelstcut(I);
	*/

	// 测试代码，枚举器
	/*
	enumerator_t en(5);

	int_set_t res;
	do
	{
		res = en.next_enumeration();
	} while (!res.empty());
	*/

	// 测试代码，枚举标签带权重的标签s-t割问题的最优解
	/*
	weighted_labelstcut_t I = generate_a_random_weighted_labelstcut_instance(
		10,	// 顶点数目
		10,	// 标签数目
		0.6,	// 两点之间有边的概率
		100	// 标签的最大权重
	);

	label_set_t opt = find_optimal_solution(I);
	*/
}


int main()
{
	test_CM();

	return 0;
}


int vertex_number(const multigraph_t& G)
{
	return G.vertex_number();
}

int edge_number(const multigraph_t& G)
{
	return G.edge_number();
}

int edge_multiplicity(const multigraph_t& G)
{
	return G.edge_multiplicity();
}


// 产生一个随机的多重图
multigraph_t generate_a_random_multigraph
(
	int n, // n -- 顶点数目
	double p, // p -- 两点之间有边的概率。
	int mu	// mu -- 两点之间尝试产生边的次数
			// 对于两个点i，j，重复mu次，每次在i和j之间以概率p产生一条边。
			// 因此，参数mu决定了边的重数。但一般地mu不一定是最大重数。
)
{
	multigraph_t G;
	if (n <= 0) return G;

	vertex_t i, j;
	int k;
	medge_set_t medge_set_empty;
		
	medge_set_empty.clear();
			
	for (i = 0; i < n; i++)
	{
		// 大括号称为“初始化列表”，这里使用了C++11的特性
		G.ve_map.insert({ i, medge_set_empty });
	}

	// 在i和j之间产生随机的边
	for (i = 0; i < n; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			// 在i和j之间尝试mu次，每次按照概率p产生一条边
			for (k = 0; k < mu; k++)
			{
				// 产生[0, 1)之间的double类型的随机数（包括0，但不包括1）
				std::uniform_real_distribution<double> dr(0, 1);
				double q = dr(dre);

				//按照概率p在i和j之间放置上边
				if (q <= p)
				{
					// e包含一个随机的id，以示i和j之间的边的区别
					medge_t e(i, j);	

					G.E.insert(e);	// 将e插入到G的边集E中
					G.ve_map[i].insert(e);	// 在顶点i的邻接表中插入e
					G.ve_map[j].insert(e);	// 在顶点j的邻接表中插入e
				}
			}
		}
	}

	return G;
}


// 按G(n, p)模型产生一个随机图
graph_t generate_a_random_graph_Gnp
(
	int n, // n -- 顶点数目
	double p // p -- 两点之间有边的概率。
)
{
	graph_t G;
	if (n <= 0) return G;
		
	vertex_t i, j;
	edge_set_t edge_set_empty;

	edge_set_empty.clear();

	// 初始化邻接表
	// 含有n个顶点，每个顶点的边集为空
	for (i = 0; i < n; i++)
	{
		// 大括号称为“初始化列表”，这里使用了C++11的特性
		G.ve_map.insert({ i, edge_set_empty });
	}

	// 产生随机的边
	for (i = 0; i < n; i++)
	{
		// 无向图的边的端点i和j是unordered的，产生随机边时，
		// 只需要尝试i之后的顶点j作为边的另一个端点即可
		for (j = i + 1; j < n; j++)
		{
			// 产生[0, 1)之间的double类型的随机数（包括0，但不包括1）
			std::uniform_real_distribution<double> dr(0, 1);
			double q = dr(dre);

			//按照概率p在i和j之间放置上边
			if (q <= p)
			{
				edge_t e(i, j);

				G.E.insert(e);	// 在边集中插入e
				G.ve_map[i].insert(e);	// 在顶点i的邻接表中插入e
				G.ve_map[j].insert(e);	// 在顶点j的邻接表中插入e
			}
		}
	}

	return G;
}


// 产生随机的度序列
int_vector_t generate_a_random_degree_sequence
(
	int n	// 顶点数目
)
{
	// 初始化整数均匀分布[1, n]
	std::uniform_int_distribution<int> di(1, n);
	// 存放顶点的度
	int_vector_t d;
	int sum = 0;

	do
	{
		d.clear();

		// 产生n个顶点的随机的度序列
		// 每个顶点的度是[1, n]之间的一个随机数
		// （度为0的情况未纳入考虑）
		for (int i = 0; i < n; i++)
			// 产生[1, n]之间的一个随机数
			d.push_back(di(dre));

		sum = 0;
		for (auto deg : d)
			sum += deg;
	} while (sum % 2 != 0);

	return d;
}


// 将d中剩余度大于0的顶点加入到c
// 将剩余度大于0的顶点和它的每个度构成的pair <i, k>加入到u
void translate_remaining_degrees
(
	int n,	// 顶点数目
	const int_vector_t& d,
	int_vector_t& c,	// 出口参数
	int_pair_vector_t& u	// 出口参数
)
{
	// 找出剩余度大于0的所有顶点
	c.clear();
	for (int i = 0; i < n; i++)
		if (d.at(i) > 0)
			c.push_back(i);

	// 将剩余度大于0的顶点和它的每个度全都加入到u
	u.clear();
	for (auto i : c)
	{
		// 现在顶点i的剩余度为d[i]
		// 将pair <i, 0>, <i, 1>, ..., <i, d[i]-1>插入到u
		for (int k = 0; k < d.at(i); k++)
		{
			u.push_back(make_pair(i, k));
		}
	}
}

// 按配置模型产生一个随机图
// 先生成一个随机的度序列d[0], d[1], ..., d[n-1]
// 再按照这个度序列，产生随机的边
// 由于过程中可能会生成自环（抛弃了）和重边，
// 最终生成的图如果计算上重边，度的序列是符合d[0], d[1], ..., d[n-1]的
// 但是，由于生成的图是简单图，最终的度序列只是接近d[0], d[1], ..., d[n-1]
graph_t generate_a_random_graph_CM
(
	int n // n -- 顶点数目
)
{
	graph_t G;
	if (n <= 0) return G;
	int_vector_t d = generate_a_random_degree_sequence(n);

	// c存放剩余度大于0的顶点
	int_vector_t c;
	// u存放顶点i和它的d[i]个度
	// 存储形式为：<i, 0>, <i, 1>, ..., <i, d[i]-1>，对每个i
	int_pair_vector_t u;

	// 将d中剩余度大于0的顶点加入到c
	// 将剩余度大于0的顶点和它的每个度构成的pair <i, k>加入到u
	translate_remaining_degrees(n, d, c, u);
	
	// 若只余下1个顶点的剩余度大于0，或没有顶点剩余度大于0，则结束
	while(c.size() >= 2)
	{
		// 初始化整数均匀分布[0, u.size()-1]
		std::uniform_int_distribution<int> di(0, u.size()-1);
		int x, y;

		// 按照度均匀分布，取两个随机的不相等的顶点
		do {
			// u.at(r)处是一个pair <i, k>，表示顶点i的第k个度
			x = u.at(di(dre)).first;
			y = u.at(di(dre)).first;
		} while (x == y);

		// 总是让x < y
		// 这样做是为了避免下面的代码将(x, y)和(y, x)都插入到E中
		// 因为edge_t并不保持i和j有序
		if(x > y) swap(x, y);

		// 在i和j之间产生一条边
		edge_t e(x, y);

		// 则将其插入到G
		// 注意：此处未做重边检查
		// 若e已经存在，则集合的insert函数不会重复插入
		G.E.insert(e);
		// 此处ve_map[x]一定要用[]运算符
		// []运算符保证了若ve_map中不存在x的map，则将x的map插入之
		G.ve_map[x].insert(e);
		G.ve_map[y].insert(e);

		// 顶点x和y的剩余度都减1
		(d.at(x))--;
		(d.at(y))--;

		// 将d中剩余度大于0的顶点加入到c
		// 将剩余度大于0的顶点和它的每个度构成的pair <i, k>加入到u
		translate_remaining_degrees(n, d, c, u);
	}

	return G;
}


string GenerateGuid()
{
	GUID guid;
	CoCreateGuid(&guid);
	char cBuffer[64] = { 0 };
	sprintf_s(cBuffer, sizeof(cBuffer), "{%08X-%04X-%04X-%02X%02X-%02X%02X%02X%02X%02X%02X}", 
		guid.Data1, guid.Data2, guid.Data3, 
		guid.Data4[0], guid.Data4[1], guid.Data4[2], guid.Data4[3], 
		guid.Data4[4], guid.Data4[5], guid.Data4[6], guid.Data4[7]);
	return string(cBuffer);
}

// 将medge_t放入set时，需要提供medge_t之间比较大小的运算符“<”的实现
bool medge_t::operator< (const medge_t &rhs) const
{
/*
	// 比较两条无向边的大小：
	// 假设两条边为e1 = (i1, j1, id1)和e2 = （i2, j2, id2)，并且假设i1 <= j1，i2 <= j2。
	// 若i1 < i2，则e1 < e2。
	// 否则，若j1 < j2，则e1 < e2。
	// 否则，若id1 < id2，则e1 < e2。
	// 否则，e1 !< e2。

	int lhs_min, rhs_min;

	lhs_min = min(i, j);
	rhs_min = min(rhs.i, rhs.j);

	if (lhs_min < rhs_min) return true;
	if (lhs_min > rhs_min) return false;

	int lhs_max, rhs_max;

	lhs_max = max(i, j);
	rhs_max = max(rhs.i, rhs.j);

	if (lhs_max < rhs_max) return true;
	if (lhs_max > rhs_max) return false;

	if (id < rhs.id) return true;

	return false;
*/

	// 先按照edge_t类型进行比较，再按照medget_t比较。
	// 若edge_t(e1) < edge_t(e2)，则e1 < e2。
	// 否则，若edge_t(e1) > edge_t(e2)，则e1 !< e2。
	// 否则（此时e1 == e2），若e1.id < e2.id，则e1 < e2。
	// 否则 e1 !< e2。
	if (edge_t(*this) < edge_t(rhs)) return true;
	if (edge_t(*this) > edge_t(rhs)) return false;

	if (id < rhs.id) return true;
	
	return false;
}


// 拷贝赋值
edge_t& edge_t::operator=(const edge_t& rhs)
{
	// 防止自赋值
	if (this == &rhs) return *this;

	i = rhs.i;
	j = rhs.j;

	return *this;
}


medge_t& medge_t::operator=(const medge_t& rhs)
{
	// 防止自赋值
	if (this == &rhs) return *this;

	edge_t::operator=(rhs);
	
	id = rhs.id;

	return *this;
}

// 图上的顶点数目
int graph_t::vertex_number() const
{
	return ve_map.size();
}

// 图上的边的数目
// 待实现
int graph_t::edge_number() const
{
	return 0;
}

// 多重图上的顶点数目
int multigraph_t::vertex_number() const
{
	return ve_map.size();
}

// 多重图上的边的数目
// 待实现
int multigraph_t::edge_number() const
{
	return 0;
}

// 多重图上的最大边重数
int multigraph_t::edge_multiplicity() const
{
	int n = vertex_number();
	int mu = 0;

	for (vertex_t i = 0; i < n; i++)
	{
		for (vertex_t j = i + 1; j < n; j++)
		{
			int m = edge_number_between_vertices(i, j);
			if (mu < m)
				mu = m;
		}
	}
		
	return mu;
}

// 返回多重图上两点之间边的数目
int multigraph_t::edge_number_between_vertices
(
	vertex_t u, 
	vertex_t v
) const
{
	int tot = 0;

	for (auto e : ve_map.at(u))
	{
		if ((e.get_i() == u && e.get_j() == v) ||
			(e.get_i() == v && e.get_j() == u))
			tot++;
	}

	return tot;
}

// 返回多重图上两点之间的边的集合
medge_set_t multigraph_t::medges_between_vertices
(
	vertex_t u, 
	vertex_t v
) const
{
	medge_set_t E_prime;

	for (auto e : ve_map.at(u))
	{
		if ((e.get_i() == u && e.get_j() == v) ||
			(e.get_i() == v && e.get_j() == u))
			E_prime.insert(e);
	}

	return E_prime;
}

// 求多重图G上从s到t的最短路，默认为无向图
// _w表示是带权重的图
// _m表示是多重图
// 返回空路径表示从s到t无路可达
mpath_t shortest_path_w_m
(
	const weighted_multigraph_t& G, 
	vertex_t s, 
	vertex_t t
)
{
	int n = vertex_number(G);

	// d存储从s到每个顶点的当前距离。
	map<vertex_t, weight_t> d;	

	// Q实现一个最小优先队列。该队列采用set实现。队列中记录的是“顶点-权重”对。
	// 之所以采用set实现优先队列，而没有采用priority_queue，
	// 是因为STL中的priority_queue是最大优先队列，而容器set中的元素是从小到大排序的。
	// 由于pair之间比较大小是按照第一部分优先的原则，因此顶点-权重对在set中
	// 实际上是按照“权重-顶点”对来记录的。
	typedef set<pair<weight_t, vertex_t>> weight_vertex_pair_set_t;
	weight_vertex_pair_set_t Q;

	typedef map<vertex_t, vertex_t> vertex_vertex_map_t;
	typedef map<vertex_t, medge_t> vertex_medge_map_t;

	// pre存储每个顶点的“前趋”（双亲）。pre实际上是一棵树。
	vertex_vertex_map_t pre;
	// branch存储顶点和其前趋（双亲）之间的边
	vertex_medge_map_t branch;
	
	// 初始化d、Q。
	for (vertex_t i = 0; i < n; i++)
	{
		// 用INT_MAX表示特别大的数M
		// 这里假设两点间的距离都小于INT_MAX
		d[i] = INT_MAX;
		Q.insert(make_pair(INT_MAX, i));
	}

	d[s] = 0;	// 将s的weight置为0
	pre[s] = s;	// s的前趋仍为s，表示s没有前趋。
	branch[s] = medge_t(s, s);	// (s, s)是一条不存在的边

	// 将顶点s的weight修改为0
	Q.erase(make_pair(INT_MAX, s));
	Q.insert(make_pair(0, s));

	/* Dijkstra算法的流程
	   从最小优先队列Q中取出队头元素u。
	   用u的权重更新V\Q中所有u的邻居的权重。
	   重复执行以上过程，直到队列为空。
	*/

	while (!Q.empty())
	{
		// 从Q中取出权重最小的顶点u
		auto p = Q.begin();
		// 调试用
		// cout << (*p).first << " " << (*p).second << endl;
		vertex_t u = (*p).second;
		Q.erase(*p);

		// 若取出权重最小的顶点后Q为空，则直接退出。
		if (Q.empty()) break;

		// 下面用u的权重去更新在Q中的u的每一个邻居v的权重
		// 首先进行检查。若d[u] == INT_MAX，则表明从s出发没有到达过u
		// 此时不需要用d[u]去更新它的在Q中的邻居的d[v]
		if (d.at(u) == INT_MAX) continue;
		
		// 下面用u的权重去更新在Q中的u的每一个邻居v的权重
		// 逐个检查u的邻居
		for (auto i = G.ve_map.at(u).begin(); i != G.ve_map.at(u).end(); ++i)
		{
			// 确定u的邻居v
			medge_t e = *i;	// (*i)的类型是medge_t
			// e的两个端点i, j中，其中一个是u。将另一个赋值给v
			vertex_t v = e.get_i();
			if (u == v) v = e.get_j();
			
			// 若v不在Q中，则不需要做任何动作
			// pos需要保存下来，下面的if中还要使用
			auto pos = Q.find(make_pair(d.at(v), v));
			if (pos == Q.end()) continue;
			
			if (d.at(u) + G.w.at(e) < d.at(v))
			{
				// 先删除掉v的老的pair
				Q.erase(pos);

				// 更新d[v]
				d.at(v) = d.at(u) + G.w.at(e);
				// 将新的pair插入到Q中
				Q.insert(make_pair(d.at(v), v));

				// 更新pre和branch
				pre[v] = u;
				branch[v] = e;
			}
		}
	}

	// p是要找到的最短路
	mpath_t p;
	// 路的起点和终点
	p.s = s;
	p.t = t;

	// 从s到t无路可达时，返回空路径
	if (d[t] == INT_MAX) return p;

	// 从s到t有路可达，则填充路的边序列，恢复出记录在pre中的最短路
	vertex_t u = t;
	while (u != s)
	{
		// 在p.e_list.begin()之前插入边branch[u]
		// branch[u]是在路上顶点u和其前趋pre[u]之间的一条边
		p.e_list.insert(p.e_list.begin(), branch.at(u));	
		u = pre.at(u);
	}
	
	return p;
}


// 求带权多重图G上路p的长度
weight_t path_length
(
	const weighted_multigraph_t& G,
	const mpath_t& p
)
{
	weight_t len = 0;

	for (auto e : p.e_list)
	{
		try
		{
			len += G.w.at(e);
		}
		catch (out_of_range)
		{
			cout << "Out of range exception. Path length may be not correct."
				<< endl;
		}
	}
	
	return len;
}


// 求带权多重图G上从s到t的距离
// 返回INT_MAX表示从s到t无路可达
weight_t dist
(
	const weighted_multigraph_t& G,
	vertex_t s,
	vertex_t t
)
{
	if (s == t) return 0;

	mpath_t p = shortest_path_w_m(G, s, t);
	if (p.e_list.empty()) return INT_MAX;

	return path_length(G, p);
}


// 产生一个随机的labelstcut_w_m_f_t实例
// 该代码尚未进行调试（2020-04-30）
labelstcut_w_m_f_t generate_a_random_labelstcut_w_m_f_instance
(
	int n,	// 顶点数目
	int q,	// 标签数目
	double p,	// 两点之间有边的概率
	int mu,	// 两个顶点之间产生边时，重复的次数
	int w_max	// 最大权重
)
{
	labelstcut_w_m_f_t I;

	// 产生随机图
	I.G = generate_a_random_multigraph(n, p, mu);

	// 默认顶点0为s，默认顶点n - 1为t
	I.s = 0;	
	I.t = n - 1;	

	// 产生q个标签
	I.A = generate_a_label_set(q);
		
	// 为G的每条边分配L中一个随机的标签

	// 尚没有好的办法产生set中随机的标签。
	// 现在采取的办法是：首先将set L复制到vector lv中
	// 由于vector中的元素都有index，就可以通过产生随机的index
	// 来产生vector中随机的标签
	label_vector_t lv = label_set_to_label_vector(I.A);
	
	// 从每个顶点开始，遍历G中所有的边
	for (vertex_t i = 0; i < n; ++i)
	{
		// 遍历顶点i所关联的边
		for (auto p = I.G.ve_map[i].cbegin(); 
			p != I.G.ve_map[i].cend(); 
			++p)
		{
			medge_t e = *p;

			// 当G是无向图时，G的整个邻接表是对称的
			// 边e = (i, j)既在i的邻接表中，也在j的邻接表中
			// 只对e产生一次随机的标签即可
			if (e.get_i() < e.get_j())
			{
				// 产生[0, lv.size()-1]之间的一个随机数
				std::uniform_int_distribution<int> di(0, lv.size() - 1);
				int k = di(dre);
				I.ell[e] = lv[k];
			}
		}
	}
	
	// 为L的每个标签分配一个随机的权重
	for (auto p = I.A.cbegin(); p != I.A.cend(); ++p)
	{
		label_t ell = *p;

		// 产生[1, w_max]之间的一个随机数
		std::uniform_int_distribution<int> di(1, w_max);
		I.w[ell] = di(dre);
	}

	return I;
}

// 产生一个随机的边上带权重的多重图
weighted_multigraph_t generate_a_random_weighted_multigraph
(
	int n,	// 顶点数目
	double p, // 两点之间有边的概率
	int mu,	// 两点之间产生边时，重复的次数
			// 对于两个点i，j，重复mu次，每次在i和j之间以概率p产生一条边。
			// 因此，参数mu决定了边的重数。但一般地mu不一定是最大重数。
	int w_max	// 最大权重
)
{
	weighted_multigraph_t G = generate_a_random_multigraph(n, p, mu);

	// 遍历邻接表的每一行，每行一个顶点
	for (auto p = G.ve_map.begin(); p != G.ve_map.end(); ++p)
	{
		// 遍历该顶点邻接的每一条边（有重边）
		for (auto q = p->second.begin(); q != p->second.end(); ++q)
		{
			medge_t e = *q;

			// 在w中查询是否已经为边e赋予了权重
			// 若e已经在w中，则跳到下一条边
			auto pos = G.w.find(e);
			if (pos != G.w.end())
				continue;

			// 产生[1, w_max]之间的一个随机数
			std::uniform_int_distribution<int> di(1, w_max);
			G.w[e] = di(dre);
		}
	}

	return G;
}


// 比较两条无向边的大小
// 比较规则如下：假设两条边为e1 = (i1, j1)和e2 = （i2, j2)
// 为简明，假设i1 <= j1，i2 <= j2（注意，边e中不保证i <= j）
// 若i1 < i2，则e1 < e2。
// 否则，若j1 < j2，则e1 < e2。
// 否则，e1 !< e2。
bool edge_t::operator< (const edge_t& rhs) const
{
	vertex_t lhs_min, rhs_min;

	lhs_min = min(i, j);
	rhs_min = min(rhs.i, rhs.j);

	if (lhs_min < rhs_min) return true;
	if (lhs_min > rhs_min) return false;

	vertex_t lhs_max, rhs_max;

	lhs_max = max(i, j);
	rhs_max = max(rhs.i, rhs.j);

	if (lhs_max < rhs_max) return true;

	return false;
}


// 比较两条无向边的大小
bool edge_t::operator> (const edge_t& rhs) const
{
	return rhs < *this;
}


//比较两条无向边是否相等
bool edge_t::operator== (const edge_t& rhs) const
{
	vertex_t lhs_min, rhs_min;

	lhs_min = min(i, j);
	rhs_min = min(rhs.i, rhs.j);

	if (lhs_min != rhs_min) return false;

	vertex_t lhs_max, rhs_max;

	lhs_max = max(i, j);
	rhs_max = max(rhs.i, rhs.j);

	if (lhs_max != rhs_max) return false;

	return true;
}


// 产生一个随机的weighted_labelstcut_t实例
// 按照ER模型中的G(n, p)模型产生随机图
weighted_labelstcut_t generate_a_random_weighted_labelstcut_instance_Gnp
(
	int n,	// 顶点数目
	int q,	// 标签数目
	double p,	// 两点之间有边的概率
	int w_max	// 标签的最大权重
)
{
	weighted_labelstcut_t I;

	// 产生随机图
	I.G = generate_a_random_graph_Gnp(n, p);

	// 默认顶点0为s，默认顶点n - 1为t
	I.s = 0;
	I.t = n - 1;

	// 为图上的边产生随机的带权重的标签
	generate_random_weighted_labels_for_edges(I, q, w_max);
		
	return I;
}


// 产生q个标签
// 标签命名法则：a, b, c, ..., z, aa, ab, ac, ..., az, ba, ...
label_set_t generate_a_label_set(int q)
{
	label_set_t L;

	// 产生q个标签
	for (int i = 0; i < q; i++)
	{
		label_t ell;
		int j = i;

		do
		{
			int k = j % 26;
			j = j / 26;
			ell = string(1, 'a' + k) + ell;
		} while (j != 0);

		L.insert(ell);
	}

	return L;
}


label_vector_t label_set_to_label_vector(const label_set_t &L)
{
	label_vector_t lv;

	for (auto p = L.cbegin(); p != L.cend(); ++p)
	{
		lv.push_back(*p);
	}

	return lv;
}

// 根据实例I中的多重图I.G生成边上带权重的多重图G
// 对I.G中的每一条边，若其标签是允许标签，则该边的权重为1
// 否则（其标签为禁止标签）该边的权重为0
weighted_multigraph_t label_weighted_multigraph_to_edge_weighted_multigraph
(
	const labelstcut_m_f_t& I
)
{
	// I.G是一个multigraph，而G是一个weighted_multigraph
	// 在weighted_multigraph_t类中定义了拷贝构造函数，支持如下赋值。
	weighted_multigraph_t G = I.G;

	for (auto e : I.G.E)
	{
		if (I.A.find(I.ell.at(e)) != I.A.end())
		{
			//若ell[e]在A中，则e的权重为1
			G.w[e] = 1;
		}
		else
		{
			G.w[e] = 0;
		}
	}

	return G;
}


// 多重图上带禁止标签的标签s-t割问题的两阶段解法
// 第一阶段，连续去掉最短路
// 返回值：第一阶段所选取的标签
// 输入参数G也发生了变化，具有选取的标签的边都被删除了
label_set_t stage_1_m_f
(
	weighted_multigraph_t& G, 
	vertex_t s, 
	vertex_t t, 
	const medge_label_map_t& ell,
	const label_set_t& A,	// 允许标签的集合
	unsigned int OPT
)
{
	int n = G.vertex_number();
	int mu = G.edge_multiplicity();	// 边的最大重数
	label_set_t sol;
	weight_t dist_bound = 
		weight_t(
			pow(n, 2.0 / 3.0) * pow(double(mu), 1.0 / 3.0) / pow(OPT, 1.0 / 3.0)
		);
	weight_t st_distance;
	
	// 当G上最短s-t路的长度<=n^2/3 mu^1/3 / OPT^1/3时
	// 连续寻找最短s-t路，并在图上删除具有p中允许标签的边
	mpath_t p = shortest_path_w_m(G, s, t);
	if (p.e_list.empty()) st_distance = INT_MAX;
	else st_distance = path_length(G, p);

	while (st_distance <= dist_bound)
	{
		// 获取L(p) \cap A中的标签
		label_set_t Lp_admissible;
		for (auto e : p.e_list)
		{
			if (find(A.begin(), A.end(), ell.at(e)) != A.end())
			{
				// 若p上的标签，同时又在A中，则插入到解中
				Lp_admissible.insert(ell.at(e));
				sol.insert(ell.at(e));
			}
		}

		// 删除G上所有具有Lp_admissible中标签的边
		remove_medges_with_specified_labels(G, ell, Lp_admissible);

		// 再计算当前最短s-t路
		p = shortest_path_w_m(G, s, t);
		if (p.e_list.empty()) st_distance = INT_MAX;
		else st_distance = path_length(G, p);
	}

	return sol;
}


// 第二阶段，对顶点分层，找一个s-t割
// 注意：一条边，若是允许边，则其长度为1。若是禁止边，则其长度为0。
// 分层是按照顶点到s的距离进行的
// 该距离并不是简单的边的数目
label_set_t stage_2_m_f
(
	weighted_multigraph_t& R,
	vertex_t s,
	vertex_t t,
	const medge_label_map_t& ell
)
{
	// 记当前图为R
	// P是对当前图R的顶点的分层
	vector<vertex_set_t> P;
	int n = R.vertex_number();
	// tau是t所在的层
	// 值为-1表示还不知到t所在的层号
	int tau = -1;	
	
	// 预先设置n个层P[0]，P[1]，...，P[n-1]
	// 并初始化每一层均为空
	for (int i = 0; i < n; i++)
	{
		vertex_set_t empty_set;
		P.push_back(empty_set);
	}

	// 遍历每一个顶点，将其放到对应的层
	for (vertex_t v = 0; v < n; v++)
	{
		int i = dist(R, s, v);
		// 距离为INT_MAX表示从s到v无路可达
		if (i == INT_MAX) continue;

		P[i].insert(v);

		// 记录t所在的层
		if (v == t) tau = i;
	}

	label_set_t sol;

	// 若tau为-1，表示没有找到t所在的层
	// 即，从s到t无路可达。此时s和t本来就是断开的，直接返回空集
	if (tau == -1) return sol;
	
	// 从s到t每相邻两层形成一个割。找出这些割中的最小者，其层号记为i_star
	int i_star = 0;
	int m_star = INT_MAX;

	for (int i = 0; i < tau; i++)
	{
		int tot = 0;

		// 计算P[i]和P[i+1]两层之间边的数目
		for (auto u : P[i])
		{
			for (auto v : P[i + 1])
			{
				tot += R.edge_number_between_vertices(u, v);
			}
		}

		// 记录当前找到的最小的层形成的割
		if (tot < m_star)
		{
			i_star = i;
			m_star = tot;
		}
	}
	
	// 将i_star层的所有边的标签都并入sol
	for (auto u : P[i_star])
	{
		for (auto v : P[i_star + 1])
		{
			label_set_t L_prime = labels(ell, R.medges_between_vertices(u, v));
			sol.insert(L_prime.begin(), L_prime.end());
		}
	}

	return sol;
}


// 求解多重图上带有禁止标签的label s-t cut问题，标签没有权重，核心算法
label_set_t core_m_f
(
	weighted_multigraph_t G, 
	vertex_t s, 
	vertex_t t, 
	const medge_label_map_t& ell, 
	const label_set_t& A,	// A是允许标签的集合
	unsigned int OPT
)
{
	// 第一阶段，连续去掉最短路
	// 注意：G也是出口参数，具有选取的标签的边都被去掉了
	label_set_t sol1 = stage_1_m_f(G, s, t, ell, A, OPT);
	
	// 第二阶段，对顶点分层，找一个s-t割
	// 注意：一条边，若是允许边，则其长度为1。若是禁止边，则其长度为0。
	// 分层是按照顶点到s的距离进行的
	// 该距离并不是简单的边的数目
	label_set_t sol2 = stage_2_m_f(G, s, t, ell);

	label_set_t sol;

	// 将sol1和sol2合并到sol
	set_union(sol1.begin(), sol1.end(), sol2.begin(), sol2.end(), 
		inserter(sol, sol.begin()));

	return sol;
}


// 求解多重图上带有禁止标签的label s-t cut问题，标签没有权重
// 返回值为求到的解（标签的集合）
// 若返回值为{"no_solution"}，则表示问题无解
label_set_t solve_labelstcut_on_multigraphs_with_forbidden_labels
(
	const labelstcut_m_f_t& I
)
{
	// 将实例I中的多重图I.G转换为边上带权重的多重图G
	// G上的边的权重仅为0或者1，分别表示相应的边为“禁止边”和“允许边”
	// 对I.G中的每一条边，若其标签是允许标签，则该边的权重为1，该边称为允许边
	// 若其标签为禁止标签，则该边的权重为0，该边称为禁止边
	weighted_multigraph_t G = label_weighted_multigraph_to_edge_weighted_multigraph(I);

	label_set_t sol;

	// 若s-t距离为0，则表示有一条路上的标签全部为forbidden标签
	// 此时返回“失败”。返回I.B表示“失败”
	if (dist(G, I.s, I.t) == 0)
	{
		sol.insert("no_solution");
		return sol;
	}

	unsigned int delta;
	// 将解初始化为I.A（选中所有允许标签）
	sol = I.A;

	// delta是对OPT的猜测
	// 对OPT的猜测值是倍增的。近似比中会有因子2
	delta = 1;
	while (delta <= I.A.size())
	{
		label_set_t sol_prime = core_m_f(G, I.s, I.t, I.ell, I.A, delta);
		// 保持sol为当前已经找到的最好的解
		if (sol_prime.size() < sol.size())
			sol = sol_prime;
		delta += delta;
	}

	// 逐个试验OPT的猜测值
	// 太耗费时间了，暂且不用
	/*
	for (delta = 1; delta <= I.A.size(); delta++)
	{
		label_set_t sol_prime = core_m_f(G, I.s, I.t, I.ell, I.A, delta);
		// 保持sol为当前已经找到的最好的解
		if (sol_prime.size() < sol.size())
			sol = sol_prime;
	}
	*/

	return sol;
}


// 将多重图G中具有L_delete中标签的边删除
// 注意G是出口参数
void remove_medges_with_specified_labels
(
	multigraph_t& G,
	const medge_label_map_t& ell,
	const label_set_t& L_delete
)
{
	// 删除G上所有具有L_delete中标签的边
	// 首先删除G.E中具有L_delete中标签的边
	remove_medges_with_specified_labels(G.E, ell, L_delete);

	// 再删除G.ve_map中具有L_delete中标签的边
	// 注意这里的auto后面是引用“&”，这是非常重要的
	// 因为在循环体中要对ve进行操作
	for (auto& ve : G.ve_map)
	{
		// ve.first是顶点
		// ve.second是边的集合（medge_set_t）
		remove_medges_with_specified_labels(ve.second, ell, L_delete);
	}
}


// 将重边的集合E中具有L中标签的边删除
// 注意E是出口参数
void remove_medges_with_specified_labels
(
	medge_set_t& E, 
	const medge_label_map_t& ell, 
	const label_set_t& L_delete
)
{
	// 将E暂存一份到E_prime
	medge_set_t E_prime = E;

	for (auto e : E_prime)
	{
		if (L_delete.find(ell.at(e)) != L_delete.end())
		{
			// 若e的标签在L_delete中，则将e从E中删除
			E.erase(e);
		}
	}
}


// 求边集所含的标签的集合
label_set_t labelstcut_m_f_t::labels(medge_set_t E)
{
	return ::labels(ell, E);
}


// 将边集E映射到标签集
label_set_t labels
(
	const medge_label_map_t& ell, 
	const medge_set_t& E
)
{
	label_set_t L_prime;

	for (medge_t e : E)
	{
		L_prime.insert(ell.at(e));
	}

	return L_prime;
}


// 求解标签s-t割问题，标签上有权重，核心算法（没有猜测步）
label_set_t core_w(labelstcut_m_f_t I)
{
	label_set_t sol;

	return sol;
}


// 返回最大的标签权重
// 若没有标签权重，则返回INT_MAX（异常情况）
weight_t weighted_labelstcut_t::max_weight()
{
	weight_set_t wset;

	for (auto i : w)
	{
		wset.insert(i.second);
	}

	if (wset.empty()) return INT_MAX;

	// max_element算法返回的是最大元的位置
	// 使用*运算符获得最大值
	return *max_element(wset.begin(), wset.end());
}


// 将ell复制w_bar_ell份，构成group_ell
label_set_t ell_to_group_ell
(
	label_t ell, 
	weight_t w_ell, 
	int q_star, 
	weight_t W
)
{
	label_set_t g_ell;

	weight_t w_bar_ell = weight_t(ceil(w_ell * double(q_star) / W));

	for (int i = 1; i <= w_bar_ell; i++)
	{
		char str[255];
		// C++编译器提示sprintf()不安全，建议使用sprintf_s()函数
		sprintf_s(str, "%d", i);
		g_ell.insert(ell + str);
	}

	return g_ell;
}


// 识别出L中的允许标签
label_set_t admissible_labels
(
	const label_set_t& L, 
	const label_weight_map_t& w, 
	weight_t W
)
{
	label_set_t A;

	// 将I.L中权重<=W的标签放入A中，是为允许标签
	for (auto l : L)
	{
		// 另一种更简洁的写法是w[l]。但在这里没有使用这种写法，
		// 原因是w[l]有副作用：若l在w中不存在，则插入之。
		// 这是map的一项特性。
		if (w.at(l) <= W)
		{
			A.insert(l);
		}
	}

	return A;
}


// 根据I、q_star和W，构造带禁止标签的多重图上的标签s-t割实例
labelstcut_m_f_t labelstcut_weighted_to_mf
(
	const weighted_labelstcut_t& I,
	int q_star,
	weight_t W
)
{
	// A是允许标签（admissible labels）的集合
	label_set_t A = admissible_labels(I.L, I.w, W);
	// B是禁止标签的集合
	label_set_t B;
	// 计算B = I.L - A
	set_difference(I.L.begin(), I.L.end(), A.begin(), A.end(), inserter(B, B.begin()));

	// 初始化生成带禁止标签的多重图上的标签s-t割实例I_tilde
	labelstcut_m_f_t I_tilde;

	// 指定I_tilde中的源和目标顶点
	I_tilde.s = I.s;
	I_tilde.t = I.t;
	
	// 定义引用变量，使用方便
	label_set_t& A_tilde = I_tilde.A;
	label_set_t& B_tilde = I_tilde.B;

	// I_tilde的禁止标签集合即是B
	B_tilde = B;

	// 下面的代码构造I_tilde的允许标签集合
	// 将A中的每个标签复制若干份，形成I_tilde的允许标签的集合A_tilde
	for (auto ell : A)
	{
		label_set_t g_ell = ell_to_group_ell(ell, I.w.at(ell), q_star, W);

		// 将A_tilde和g_ell合并到A_tilde
		label_set_t t;
		set_union(A_tilde.begin(), A_tilde.end(), g_ell.begin(), g_ell.end(),
			inserter(t, t.begin()));
		A_tilde = t;
	}

	char str[255];
	// 定义引用变量，使用方便
	multigraph_t& G_tilde = I_tilde.G;

	// 初始化G_tilde的ve_map
	int n = I.G.vertex_number();

	for (int i = 0; i < n; i++)
	{
		medge_set_t empty_set;
		G_tilde.ve_map[i] = empty_set;
	}

	// 下面根据G生成G_tilde，以及I_tilde.ell
	// 遍历G中的每一条边e，若其标签是允许标签，则在G_tilde中生成w_bar_ell数目的重边
	// 每条重边上的标签都是新实例中的允许标签
	// 若其标签是禁止标签，则直接将其复制到G_tilde（不生成重边）
	// 在新的实例中，其上的标签仍为禁止标签
	for (auto e : I.G.E)
	{
		if (find(A.begin(), A.end(), I.ell.at(e)) != A.end())
		{
			// 若ell[e]在A中，则ell[e]是允许标签
			int w_bar_ell = int(ceil(I.w.at(I.ell.at(e)) * double(q_star) / W));
			for (int i = 1; i <= w_bar_ell; i++)
			{
				// 生成e的拷贝e_copy。e_copy的id随机自动生成
				medge_t e_copy(e.get_i(), e.get_j());
				// 将e_copy放入G_tilde的边集E中
				G_tilde.E.insert(e_copy);
				// 将e_copy分别放入其两个顶点的邻接表中
				G_tilde.ve_map[e_copy.get_i()].insert(e_copy);
				G_tilde.ve_map[e_copy.get_j()].insert(e_copy);
								
				// 给e_copy分配标签ell_i
				sprintf_s(str, "%d", i);
				label_t ell_i = I.ell.at(e) + str;
				// ell_i应该在A_tilde中。检查确认
				if (find(A_tilde.begin(), A_tilde.end(), ell_i) == A_tilde.end())
				{
					cout << "Fatal error: " << endl;
					exit(1);
				}
				I_tilde.ell[e_copy] = ell_i;
			}
		}
		else
		{
			// 若ell[e]不在A中，则ell[e]是禁止标签
			// 仅将e的一份拷贝放入G_tilde
			medge_t e_copy(e.get_i(), e.get_j());
			G_tilde.E.insert(e_copy);
			// 将e_copy分别放入其两个顶点的邻接表中
			G_tilde.ve_map[e_copy.get_i()].insert(e_copy);
			G_tilde.ve_map[e_copy.get_j()].insert(e_copy);
			// 给e_copy分配标签
			I_tilde.ell[e_copy] = I.ell.at(e);
		}
	}

	return I_tilde;
}


// 判断S是否为A的子集
// 需要进行测试（2020-08-20）
bool is_label_subset
(
	const label_set_t& S, 
	const label_set_t& A
)
{
	label_set_t r1, r2;

	// 计算r1 = A - S
	set_difference(A.begin(), A.end(), S.begin(), S.end(), inserter(r1, r1.begin()));
	// 计算r2 = A /_\ S
	set_symmetric_difference(A.begin(), A.end(), S.begin(), S.end(), 
		inserter(r2, r2.begin()));

	return r1 == r2;
}

// 将group labels的集合转换为labels的集合
label_set_t group_labels_to_labels
(
	const label_set_t& A_tilde_prime,
	const label_set_t& L,
	const label_weight_map_t& w,
	int q_star,
	weight_t W
)
{
	label_set_t L_double_prime;

	for (auto ell : L)
	{
		if (w.at(ell) <= W)
		{
			label_set_t g_ell = ell_to_group_ell(ell, w.at(ell), q_star, W);
			// 判断g_ell是否为A_tilde_prime的子集
			if (is_label_subset(g_ell, A_tilde_prime))
			{
				L_double_prime.insert(ell);
			}
		}
	}

	return L_double_prime;
}

// 求解标签带权重的label s-t cut问题
label_set_t solve_weighted_labelstcut(const weighted_labelstcut_t& I)
{
	// q_star为猜测的最优解所使用的标签的数目
	int q_star;
	
	// wset用于存储所有标签的不同权重，是猜测W的取值集合
	weight_set_t wset;

	// w是一个label_weight_map_t，其每一个元素是一个pair
	// pair的第一部分是标签，pair的第二部分是权重
	for (auto i : I.w)
		wset.insert(i.second);

	// L_prime是找到的解
	label_set_t L_prime = I.L;

	// 尝试q_star的每一个可能的值，是为猜测
	for (q_star = 1; size_t(q_star) <= I.L.size(); q_star++)
	{
		// W为猜测的最优解中的最大标签权重
		// 尝试W的每一个可能的值，是为猜测
		for (auto W : wset)
		{
			// 根据I、q_star和W，构造带禁止标签的多重图上的标签s-t割实例
			labelstcut_m_f_t I_tilde = labelstcut_weighted_to_mf(I, q_star, W);
			// 求解带禁止标签的多重图上的标签s-t割实例
			label_set_t A_tilde_prime = 
				solve_labelstcut_on_multigraphs_with_forbidden_labels(I_tilde);
			// 若没有找到解，则尝试下一个W
			if (A_tilde_prime.find("no_solution") != A_tilde_prime.end()) 
				continue;

			// 将求解结果转换为标签的集合
			label_set_t L_double_prime = 
				group_labels_to_labels(A_tilde_prime, I.L, I.w, q_star, W);
			// 若找到更好的解，则用其替换当前的解
			if (L_double_prime.size() < L_prime.size())
			{
				L_prime = L_double_prime;
			}
		}
	}
	
	return L_prime;
}


void enumerator_t::initialize(int n) 
{ 
	this->n = n; 
	current_enum.clear();
}

int_set_t enumerator_t::next_enumeration()
{
	if (n <= 0)
	{
		current_enum.clear();
		return current_enum;
	}

	// 判断是否为第一次枚举
	if (current_enum.size() == 0)
	{
		current_enum.insert(0);
		return current_enum;
	}

	// 根据当前枚举current_enum生成下一个枚举
	// 先将current_enum放到vector中
	int_vector_t c;

	for (auto e : current_enum)
	{
		c.push_back(e);
	}

	// 定位到c的最后一个元素
	int len = c.size();
	int i = len - 1;
	bool found = false;

	while (!found)
	{
		// 判断是否可以将c[i]加1，并且c[i]之后的元素的值依次递增
		if (c[i] + 1 + (len - 1 - i) <= n - 1)
		{
			// 将c[i]加1，c[i]之后的元素的值依次递增
			c[i]++;

			for (++i; i < len; i++)
				c[i] = c[i - 1] + 1;
			found = true;
			continue;
		}
		else
		{
			// 无法递增，则回退i
			if (i > 0)
			{
				// 可以回退
				// 然后返回循环顶部，继续尝试增长
				i--;
				continue;
			}
			else
			{
				// 若i已经等于0，则无法回退
				if (len < n)
				{
					// 若长度 < n，则让序列长度加1
					// 再去枚举序列长度加1之后的第一个序列
					c.push_back(0);
					len++;
					// 为c[i]加1做好准备
					c[i] = -1;
					continue;
				}
				else
					// 若序列长度已经等于n，则已经枚举到了最后一个序列
					break;
			}
		}
	}

	current_enum.clear();

	if(found)
	{ 
		for (auto e : c)
		{
			current_enum.insert(e);
		}
	}
	
	return current_enum;
}


// 使用枚举的方法，求标签带权重的标签st割问题的最优解
label_set_t find_optimal_solution
(
	const weighted_labelstcut_t& I
)
{
	label_set_t sol;
	if (!connected(I.G, I.s, I.t))
		return sol;

	enumerator_t en(I.L.size());
	label_vector_t v = label_set_to_label_vector(I.L);
	sol = I.L;

	//生成第一个枚举（下标的枚举）
	int_set_t seq = en.next_enumeration();
	while (!seq.empty())
	{
		// 按照seq指定的下标，从v中取出对应标签，构成sol_prime
		label_set_t sol_prime;
		for (auto i : seq)
			sol_prime.insert(v[i]);

		if (weight(sol_prime, I.w) < weight(sol, I.w))
		{
			// 仅当sol_prime的权重小于sol时，才尝试其是不是解
			graph_t G_prime = I.G;
			remove_edges_with_specified_labels(G_prime, I.ell, sol_prime);
			if (!connected(G_prime, I.s, I.t))
			{
				// sol_prime是解，且权重小于sol，更新sol
				sol = sol_prime;
			}
		}
		
		// 取下一个枚举（下标的枚举）
		seq = en.next_enumeration();
	}

	return sol;
}


// 将图G中具有L_delete中标签的边删除
// 注意G是出口参数
void remove_edges_with_specified_labels
(
	graph_t& G,
	const edge_label_map_t& ell,
	const label_set_t& L_delete
)
{
	// 删除G上所有具有L_delete中标签的边
	// 首先删除G.E中具有L_delete中标签的边
	remove_edges_with_specified_labels(G.E, ell, L_delete);

	// 再删除G.ve_map中具有L_delete中标签的边
	// 注意这里的auto后面是引用“&”，这是非常重要的
	// 因为在循环体中要对ve进行操作
	for (auto& ve : G.ve_map)
	{
		// ve.first是顶点
		// ve.second是边的集合（medge_set_t）
		remove_edges_with_specified_labels(ve.second, ell, L_delete);
	}
}


// 将边集E中具有L中标签的边删除
// 注意E是出口参数
void remove_edges_with_specified_labels
(
	edge_set_t& E,
	const edge_label_map_t& ell,
	const label_set_t& L_delete
)
{
	// 将E暂存一份到E_prime
	edge_set_t E_prime = E;

	// 检查E_prime的每一条边，若它的标签在L_delete中，则将该边从E中删除
	for (auto e : E_prime)
	{
		if (L_delete.find(ell.at(e)) != L_delete.end())
		{
			// 若e的标签在L_delete中，则将e从E中删除
			E.erase(e);
		}
	}
}


// 判断在G中s和t是否连通
bool connected
(
	const graph_t& G,
	vertex_t s,
	vertex_t t
)	
{
	return BFS_conditional(G, s, t);
}


// 从s开始BFS，遇到t则结束
// 返回值：BFS过程中遇到了t则返回真，否则返回假
bool BFS_conditional
(
	const graph_t& G, 
	vertex_t s,
	vertex_t t
)
{
	if (s == t) return true;

	int n = G.vertex_number();
	bool_vector_t visited(n);
	queue<vertex_t> Q;

	Q.push(s);

	for (int i = 0; i < n; i++)
		visited[i] = false;

	while (!Q.empty())
	{
		vertex_t u = Q.front();
		Q.pop();

		if (!visited[u])
		{
			visited[u] = true;

			for (auto e : G.ve_map.at(u))
			{
				int v = e.get_i();
				if (v == u) v = e.get_j();

				if (v == t) return true;

				Q.push(v);
			}
		}
	}

	return false;
}


// 计算标签集的权重
weight_t weight
(
	const label_set_t& L,
	const label_weight_map_t& w
)
{
	weight_t w_tot = 0;

	for (auto ell : L)
	{
		w_tot += w.at(ell);
	}

	return w_tot;
}


// 以G(n, p)模型产生随机实例
// 测试标签带权重的标签s-t割问题的近似解法的性能
void test_Gnp()
{
	weighted_labelstcut_t I;
	label_set_t opt;
	label_set_t sol;
	weight_t OPT = 0;
	weight_t SOL = 0;

	// 同一个参数，产生10个随机实例
	// 求OPT和SOL的平均值
	for (int i = 0; i < 10; i++)
	{
		// 以G(n, p)模型产生随机实例
		I = generate_a_random_weighted_labelstcut_instance_Gnp(
			6,	// 顶点数目
			10,	// 标签数目
			0.6,	// 两点之间有边的概率
			100	// 标签的最大权重
		);

		opt = find_optimal_solution(I);
		sol = solve_weighted_labelstcut(I);
		OPT += weight(opt, I.w);
		SOL += weight(sol, I.w);
		cout << i << "\t";
		cout << "OPT = " << weight(opt, I.w) << "\t";
		cout << "SOL = " << weight(sol, I.w) << endl;
	}
	
	cout << "Average OPT = " << OPT/10 << "\t";
	cout << "Average SOL = " << SOL/10 << endl;
}


// 以配置模型产生随机实例
// 测试标签带权重的标签s-t割问题的近似解法的性能
void test_CM()
{
	weighted_labelstcut_t I;
	label_set_t opt;
	label_set_t sol;
	weight_t OPT = 0;
	weight_t SOL = 0;

	// 同一个参数，产生10个随机实例
	// 求OPT和SOL的平均值
	for (int i = 0; i < 10; i++)
	{
		// 以G(n, p)模型产生随机实例
		I = generate_a_random_weighted_labelstcut_instance_CM(
			12,	// 顶点数目
			16,	// 标签数目
			100	// 标签的最大权重
		);

		opt = find_optimal_solution(I);
		sol = solve_weighted_labelstcut(I);
		OPT += weight(opt, I.w);
		SOL += weight(sol, I.w);
		cout << i << "\t";
		cout << "OPT = " << weight(opt, I.w) << "\t";
		cout << "SOL = " << weight(sol, I.w) << endl;
	}

	cout << "Average OPT = " << OPT / 10 << "\t";
	cout << "Average SOL = " << SOL / 10 << "\t";
	cout << "ratio = " << double(SOL) / OPT << endl;

}


// 产生一个随机的weighted_labelstcut_t实例
// 按照配置模型产生随机图
weighted_labelstcut_t generate_a_random_weighted_labelstcut_instance_CM
(
	int n,	// 顶点数目
	int q,	// 标签数目
	int w_max	// 标签的最大权重
)
{
	weighted_labelstcut_t I;

	// 产生随机图
	I.G = generate_a_random_graph_CM(n);

	// 默认顶点0为s，默认顶点n - 1为t
	I.s = 0;
	I.t = n - 1;

	// 为边产生随机的带权重的标签
	generate_random_weighted_labels_for_edges(I, q, w_max);

	return I;
}

// 为L中的每个标签产生一个随机的权重，权重取值为[1, w_max]
label_weight_map_t generate_label_weights
(
	const label_set_t& L, 
	weight_t w_max
)
{
	label_weight_map_t w;

	// 为L的每个标签分配一个随机的权重
	for (auto ell:L)
	{
		// 产生[1, w_max]之间的一个随机数
		std::uniform_int_distribution<int> di(1, w_max);
		w[ell] = di(dre);
	}

	return w;
}


// 为图上的边产生随机的带权重的标签
// 分为如下步骤：
// （1）产生q个标签，并为标签产生随机的权重
// （2）为每条边分配一个随机的标签
void generate_random_weighted_labels_for_edges
(
	weighted_labelstcut_t& I,
	int q,
	weight_t w_max
)
{
	int n = I.G.vertex_number();

	// 产生q个标签
	I.L = generate_a_label_set(q);

	// 为每个标签产生一个随机的权重，权重取值为[1, w_max]
	I.w = generate_label_weights(I.L, w_max);

	// 为G的每条边分配L中一个随机的标签
	//
	// 尚没有好的办法产生set中随机的标签。
	// 现在采取的办法是：首先将set L复制到vector lv中
	// 由于vector中的元素都有index，就可以通过产生随机的index
	// 来产生vector中随机的标签
	label_vector_t lv = label_set_to_label_vector(I.L);

	// 从每个顶点开始，遍历G中所有的边
	for (vertex_t i = 0; i < n; ++i)
	{
		// 遍历顶点i所关联的边
		for (auto p = I.G.ve_map[i].cbegin();
			p != I.G.ve_map[i].cend();
			++p)
		{
			edge_t e = *p;

			// 由于G是无向图，G的整个邻接表是对称的
			// 边e = (i, j)既在i的邻接表中，也在j的邻接表中
			// 只对e产生一次随机的标签即可
			if (e.get_i() < e.get_j())
			{
				// 产生[0, lv.size()-1]之间的一个随机数
				std::uniform_int_distribution<int> di(0, lv.size() - 1);
				int k = di(dre);
				I.ell[e] = lv[k];
			}
		}
	}

	// 由于边上的标签是随机选取的，并不能保证L中的标签都出现在了E中
	// 下面的代码“删除”L和w中没有使用的标签
	label_set_t L_prime;

	for (auto e : I.G.E)
	{
		L_prime.insert(I.ell.at(e));
	}

	// U是没有使用的标签的集合
	label_set_t U;

	// 计算U = L - L_prime
	set_difference(I.L.begin(), I.L.end(), L_prime.begin(), L_prime.end(),
		inserter(U, U.begin()));

	for (auto ell : U)
	{
		// 删除map中的ell项
		I.w.erase(ell);
	}

	// 仅保留使用过的标签
	I.L = L_prime;
}

// EOF: lc.cpp