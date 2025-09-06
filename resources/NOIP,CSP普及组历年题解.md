# 2003

## T1 乒乓球-模拟
```cpp
/*
11分制：
1）分差大于2
2）某一方比分大于等于11

21分制度：
1）分差大于2
2）某一方比分大于等于21

E不一定只出现在末尾，有可能有这种情况：WLEWWWWW，E之后内容要忽略。
*/
#include<iostream>
#include<string>
#define io(x); freopen(x".in", "r", stdin), freopen(x".out", "w", stdout);
using namespace std;

string s = "", k;
int p1, p2;

void solve(int type)
{
	p1 = 0, p2 = 0;
	for (int i = 0; s[i] != 'E'; i++)
	{
		if (s[i] == 'W') p1++;
		if (s[i] == 'L') p2++;
		if ((p1 >= type || p2 >= type) && abs(p1-p2) > 1)
		{
			cout << p1 << ":" << p2 << endl;
			p1 = 0;
			p2 = 0;
		}
	}
	cout << p1 << ":" << p2 << endl;
}

int main()
{
	io("pingpong");
	while(cin >> k) s += k;
	solve(11);
	cout << endl;
	solve(21);
	return 0;
}
```
## T2 数字游戏-区间dp，前缀和
```cpp
#include<bits/stdc++.h>
using namespace std;
const int N = 110, M = 10, INF = 0x3f3f3f3f;

int n, m;
int w[N], s[N];
int f[N][N][M], g[N][N][M]; // f/g[i][j][k]表示将i~j区间分割成k部分的最小/大和

int get_mod(int x) // 处理负数取模
{
    return (x % 10 + 10) % 10;
}

int main()
{
    freopen("game.in","r",stdin);
    freopen("game.out","w",stdout);
    cin >> n >> m;
    for (int i = 1; i <= n; i ++ )
    {
        cin >> w[i];
        w[i + n] = w[i]; // 破环成链
    }
    for (int i = 1; i <= n * 2; i ++ ) s[i] = s[i - 1] + w[i]; // 求前缀和

    memset(f, INF, sizeof f); // 初始化极大值
    memset(g, -INF, sizeof g); // 初始化极小值
	// 区间dp模板，一层枚举区间长度，二层枚举区间起点，三层枚举分割点/数
    for (int len = 1; len <= n; len ++ )
        for (int l = 1; l + len - 1 <= n * 2; l ++ )
        {
            int r = l + len - 1;
            f[l][r][1] = g[l][r][1] = get_mod(s[r] - s[l - 1]);
            for (int k = 2; k <= m; k ++ )
                for (int j = l + k - 2; j < r; j ++ )
                {
                    f[l][r][k] = min(f[l][r][k], f[l][j][k - 1] * get_mod(s[r] - s[j]));
                    g[l][r][k] = max(g[l][r][k], g[l][j][k - 1] * get_mod(s[r] - s[j]));
                }
        }

    int maxv = -INF, minv = INF;
    for (int i = 1; i <= n; i ++ )
    {
        maxv = max(maxv, g[i][i + n - 1][m]);
        minv = min(minv, f[i][i + n - 1][m]);
    }

    cout << minv << endl;
    cout << maxv << endl;

    return 0;
}
```

# 2004
## T1 不高兴的津津-模拟
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	int maxday = -1, day = 0;
	for (int i = 1; i <= 7; i++)
	{
		int a, b;
		cin >> a >> b;
		if (a + b > maxday && a + b > 8) // 找>8的最小的i
			maxday = a + b, day = i;
	}
	cout << day;
	
	return 0;
}
```

## T2 花生采摘-模拟，结构体排序
```cpp
#include<bits/stdc++.h>
using namespace std;

struct Node
{
	int x, y; // 坐标
	int val; // 花生个数
	int t;  // 采完当前花生需要的时间
};
Node p[450]; // 存所有花生的信息
int a[25][25];

bool cmp(Node a, Node b)
{
	return a.val > b.val;
}

int main()
{
	freopen("hs.in", "r", stdin);
	freopen("hs.out", "w", stdout);
	int m, n, k, idx = 1;
	cin >> m >> n >> k;
	for (int i = 1; i <= m; i++) // 数据读入
		for (int j = 1; j <= n; j++)
		{
			cin >> a[i][j];
			if (a[i][j] != 0)
			{
				p[idx].x = i;
				p[idx].y = j;
				p[idx].val = a[i][j];
				idx ++;
			}
		}
	
	// 按花生由大到小排序
	sort(p+1, p+1+idx, cmp);
	// 模拟采摘过程，由大到小遍历，判断采完下一个花生可以回到路边，就去采
	int ans = 0;
	for (int i = 1; i <= idx; i++)
	{
		if (i == 1) p[i].t = p[i].x + 1; // 马路跳到田里需要额外1单位时间
		else p[i].t = p[i-1].t + abs(p[i-1].x-p[i].x) + abs(p[i-1].y-p[i].y) + 1;
		if (k - p[i].t >= p[i].x) ans += p[i].val;
	}
	cout << ans;
	return 0;
}
```

# 2005

## T1 陶陶摘苹果-模拟
```cpp
#include<bits/stdc++.h>  
using namespace std;  
int a[15];  
int main()  
{  
    int cnt = 0;  
    for (int i = 1; i <= 10; i++) cin >> a[i];  
    int h;  
    cin >> h;  
    for (int i = 1; i <= 10; i++)  
    {  
        if (h + 30 >= a[i])  
            cnt++;  
    }  
    cout << cnt;  
    return 0;  
}
```

## T2 校门外的树-数组计数
```cpp
#include<iostream>  
using namespace std;  

const int N = 1e4+10;  
int st[N];  

int main()  
{  
	int L, M;  
	cin >> L >> M;  
	while (M--)  
	{  
		int l, r; // [l,r]挖走  
		cin >> l >> r;  
		for (int i = l; i <= r; i++) st[i] = 1;  
	}  
	int cnt = 0;  
	for (int i = 0; i <= L; i++)  
		if (st[i] == 0) cnt ++;  
	cout << cnt;  
	return 0;  
}
```

# 2006
## T1 明明的随机数-桶排序
```cpp
#include<bits/stdc++.h>
using namespace std;

set<int> s;
int main()
{
	int n;
	cin >> n;
	for (int i = 1; i <= n; i++)
	{
		int x;
		cin >> x;
		s.insert(x);
	}
	cout << s.size() << endl;
	for (auto it : s)
		cout << it << " ";
	
	return 0;
}
```

## T2 开心的金明-01背包问题
```cpp
/*
总钱n相当于背包容量
每件物品的加个相当于体积
每件物品的加个乘以重要度相当于价值
*/
#include<bits/stdc++.h>
using namespace std;
const int N = 30010;
int n, m;
int dp[N];

int main()
{
	cin >> m >> n;
	for (int i = 0; i < n; i ++ )
	{
		int v, w;
		cin >> v >> w;
		w *= v;
		for (int j = m; j >= v; j -- )
			dp[j] = max(dp[j], dp[j - v] + w);
	}
	cout << dp[m] << endl;
	
	return 0;
}
```

# 2007
## T1 奖学金-结构体排序
```cpp
#include <iostream>
#include <algorithm> 

using namespace std;
struct grade // 定义结构体 
{
	int num;
	int cn, math, es;
	int sum;
};

grade stu[310]; // 生成一个结构体对象 

bool compare(grade a, grade b) // 定义比较函数 
{	// 总分大的在前，总分相同语文大的在前，总分语文都相同学号小的在前 
	if (a.sum != b.sum) return a.sum > b.sum; // 总分大的在前 
	else
	{
		if (a.cn != b.cn) 
			return a.cn > b.cn; // 语文大的在前 
		else return a.num < b.num; // 编号小的在前 
	} 
}

int main()
{	
	int n;
	cin >> n; // 输入
	for (int i = 1; i <= n; i++)
	{
		cin >> stu[i].cn >> stu[i].math >> stu[i].es;
		stu[i].num = i;
		stu[i].sum = stu[i].cn + stu[i].math + stu[i].es;	
	} 
	sort(stu+1, stu+n+1, compare);
	for(int i = 1; i <= 5; i++)
		cout << stu[i].num << " " << stu[i].sum << endl;
	
	return 0;
}
```

## T2 纪念品分组-贪心，双指针
```cpp
/*
要使得组数尽可能少，应该让每组的价值之和尽可能大，由此有如下思路：
1.将所有物品按价值升序排序
2.最小和最大匹配，判定是否可以分组，可以则考虑次小和次大匹配，不可以则最大的只能单独一组。
*/

#include<bits/stdc++.h>
using namespace std;
const int N = 30010;
int n, m;
int w[N];

int main()
{
	freopen("group.in", "r", stdin);
	freopen("group.out", "w", stdout);
	cin >> m >> n;
	for (int i = 0; i < n; i ++ ) cin >> w[i];
	
	sort(w, w + n);
	
	int res = 0;
	int i = 0, j = n-1; // 碰撞指针指向头尾
	while (i < j)
	{
		if (w[i] + w[j]  <= m) // 最小最大划分一组
			res++, i++, j--;
		else
			res++, j--; // 最大单独一组
	}
	if (i == j) res++; // 若有剩下，剩下的单独一组
	cout << res << endl;
	return 0;
}
```

# 2008
## T1 ISBN号-模拟，字符串处理
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	freopen("isbn.in", "r", stdin);
	freopen("isbn.out", "w", stdout);
	string str;
	cin >> str;
	int sum = 0;
	for (int i = 0, j = 1; i + 1 < str.size(); i ++ )
		if (str[i] != '-')
		{
			sum += (str[i] - '0') * j;
			j ++ ;
		}
	sum %= 11;
	char c = 'X';
	if (sum < 10) c = '0' + sum;
	if (c == str.back()) cout << "Right";
	else
	{
		str.back() = c;
		cout << str << endl;
	}
	return 0;
}
```

## T2 排座椅-贪心，排序
```cpp
/*  
这道题目的核心，是需要发现如下性质：  
  
不同行、列之间是完全独立的。即不管将哪行、哪列切开，对其余的行列都是没有任何影响的。  
  
因此可以分别考虑行和列。  
对于行来说，问题变成：  
  
去掉哪 K行，可以使得最后剩下的行间的边数最少。这里去掉边数最多的 K  
行一定是最优的。否则可以将选出的行替换成边数最多的 K 行，且结果不会变差。  
  
*/  
  
#include<bits/stdc++.h>  
using namespace std;  
typedef pair<int, int> PII;  
const int N = 1010;  
  
int n, m, L, K, D;  
PII row[N], col[N];  
int q[N];  
  
int main()  
{  
    cin >> n >> m >> K >> L >> D;  
        for (int i = 1; i <= n; i ++ ) row[i].second = i;  
    for (int i = 1; i <= m; i ++ ) col[i].second = i;  
        while (D -- )  
    {  
        int x1, y1, x2, y2;  
        cin >> x1 >> y1 >> x2 >> y2;  
        if (abs(x1 - x2) == 1) row[min(x1, x2)].first ++ ;  
        else col[min(y1, y2)].first ++ ;  
    }  
        sort(row + 1, row + n + 1);  
    sort(col + 1, col + m + 1);  
        int cnt = 0;  
    for (int i = n; i > n - K; i -- ) q[cnt ++ ] = row[i].second;  
    sort(q, q + cnt);  
    for (int i = 0; i < cnt; i ++ ) cout << q[i] << " ";  
    cout << endl;  
        cnt = 0;  
    for (int i = m; i > m - L; i -- ) q[cnt ++ ] = col[i].second;  
    sort(q, q + cnt);  
    for (int i = 0; i < cnt; i ++ ) cout << q[i] << " ";  
    cout << endl;  
        return 0;  
}
```

# 2009
## T1 多项式输出-模拟
```cpp
#include <bits/stdc++.h>
#define io(x); freopen(x".in", "r", stdin), freopen(x".out", "w", stdout);
using namespace std;

int main()
{
	io("poly");
	int n;
	cin >> n;
	bool is_first = true;
	for (int i = n; i >= 0; i -- )
	{
		int a;
		cin >> a;
		if (!a) continue; // 系数0，不输出
		if (!is_first && a > 0) printf("+"); // 非首位正数，输出+
		else if (a < 0) printf("-"); // 负数输出-
		if (abs(a) != 1 || !i) printf("%d", abs(a)); // 系数非1，且幂非0，输出系数绝对值
		if (i) printf("x"); // 幂指数非0，x^存在
		if (i > 1) printf("^%d", i); // 幂指数 > 1，保留指数
		
		is_first = false;
	}
	return 0;
}
```

## T2 分数线划定-模拟，排序
按分数高，序号小的排序，求出排名在$1.5\times m$的分数作为分数线，找出分数大于分数线的同学
```cpp
#include<bits/stdc++.h>
using namespace std;

const int N = 5010;
struct Node
{
	int id, scoce;
}p[N];

bool cmp(Node a, Node b)
{
	if (a.scoce == b.scoce) return a.id < b.id;
	return a.scoce > b.scoce;
}

int main()
{
	int n, m;
	cin >> n >> m;
	
	for (int i = 1; i <= n; i++)
		cin >> p[i].id >> p[i].scoce;
	sort(p+1, p+1+n, cmp);
	
	int line = p[int(m*1.5)].scoce; // 分数线
	int cnt = 0;
	for (int i = 1; i <= n; i++)
		if (p[i].scoce >= line)
			cnt++;
	cout << line << " " << cnt << endl;
	for (int i = 1; i <= n; i++)
		if (p[i].scoce >= line)
			cout << p[i].id << " " << p[i].scoce << endl;
	return 0;
}
```

# 2010
## T1 数字统计-模拟，数位分离
```cpp
#include<bits/stdc++.h>
using namespace std;

int get2(int x)
{
	int res = 0;
	while (x)
	{
		if (x % 10 == 2) res++;
		x /= 10;
	}
	return res;
}

int main()
{
	freopen("tj.in", "r", stdin);
	freopen("tj.out", "w", stdout);
	int l, r;
	cin >> l >> r;
	int cnt = 0;
	for (int i = l; i <= r; i++) cnt += get2(i);
	cout << cnt;
	
	return 0;
}
```

## T2 接水问题-贪心，模拟
模拟接水过程，每次找出当前接水量最小的水龙头，将第i个同学安排过去接水，水龙头m<=100，直接循环遍历查找即可，也可用优先队列进行优化。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
const int N = 10010, M = 110;  
int n, m;  
int w[N];  
int q[M];  
  
int main()  
{  
    cin >> n >> m;  
    for (int i = 0; i < n; i ++ ) cin >> w[i];  
    for (int i = 0; i < n; i ++ )  
    {  
        int t = 0;  
        for (int j = 0; j < m; j ++ ) // 求最小龙头  
            if (q[j] < q[t])  
                t = j;  
        q[t] += w[i]; // 第i个同学续上  
    }  
    int res = 0;  
    for (int i = 0; i < m; i ++ )         res = max(res, q[i]);  
        cout << res << endl;  
        return 0;  
}
```

# 2011
## T1 数字反转-模拟
（1）数位分离
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	int n, rev_n = 0;
	cin >> n;
	while (n)
	{
		rev_n = rev_n * 10 + n % 10;
		n /= 10;
	}
	cout << rev_n;
	return 0;
}
```
（2）字符串
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	string s;
	cin >> s;
	bool f = 0;
	if (s[0] == '-')
	{
		f = 1;
		s = s.substr(1);
	}
	while (s.size() > 1 && s.back() == '0') s.pop_back();
	reverse(s.begin(), s.end());
	if (f) s = '-' + s;
	cout << s;
	return 0;
}
```

## T2 统计单词数-模拟，字符串
读入待查找单词t和查找目标s，双指针找到每个单词k和t进行对比，记录首个位置和数量即可。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
void tolow(string &s) // 转小写  
{  
    for (auto &it:s)  
        it = tolower(it);  // tolower为内置函数，可以将单字符转小写
}  
  
int main()  
{  
    freopen("cal.in", "r", stdin);  
    freopen("cal.out", "w", stdout);  
    string t, s;  
    getline(cin, t); // 都需要用getline，cin>>t会导致无法完整读入s  
    getline(cin, s);  
    tolow(t), tolow(s);  
    int cur = -1, cnt = 0;  
    for (int i = 0; i < s.size(); i++)  
    {  
        int j = i;  
        while (j < s.size() && s[j] != ' ') j++;  
        string k = s.substr(i, j-i); // 句中单词k
        if (t == k)  
        {  
            if (cur == -1) cur = i;  
            cnt++;  
        }  
        i = j;  
    }  
    if (!cnt) cout << -1;  
    else cout << cnt << " " << cur;  
    return 0;  
}
```

# 2012
## T1 质因数分解-数学
题目明确说了n是两个不同质数的乘积，从2开始枚举的第一个因子一定是较小的那个质因子，不需要进行质数判定。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    int n;  
    cin >> n;  
    for (int i = 2; i <= n / i; i++)  
    {  
        if (n % i == 0)  
        {  
            cout << n / i;  
            break;  
        }  
    }  
    return 0;  
}
```

## T2 寻宝-模拟
本题直接模拟会超时，因为房间指示牌数字最大1e6，相当于在同一楼层一直转圈浪费时间，仔细观察发现每层房间数最多100，只需要统计出当前层有楼梯的房间数量`l`，用指示牌数字`%l`即可完成优化。

时间复杂度：
优化前：O(nx)
优化后：O(nm)
```cpp
#include<bits/stdc++.h>  
using namespace std;  
const int N=10010,M=110;  
int s[N][M],x[N][M];//一个存有没有楼梯，一个存指示牌得数  
int main(){  
    int n,m,k;  
    scanf("%d%d",&n,&m);  
    for(int i=0;i<n;i++)for(int j=0;j<m;j++)scanf("%d%d",&s[i][j],&x[i][j]);//双层循环输入  
    scanf("%d",&k);//输入第一个进的房间  
    int r=0;  
    for(int i=0;i<n;i++){  
        int t=x[i][k];  
        int l=0;//存本层楼梯数量  
        for(int j=0;j<m;j++)l+=s[i][j];//计算本层楼梯数量  
        t%=l;//取余优化  
        if(t==0)t=l;//特判  
        r=(r+x[i][k])%20123;//计算秘钥  
        for(int j=k;;j=(j+1)%m){//计算第s[i][k]个有楼梯的房间  
            if(s[i][j]){//当层是否有楼梯  
                t--;  
                if(t==0){//t是否为零  
                    k=j;//上楼  
                    break;//进入下一个循环  
                }  
            }  
        }  
    }  
    printf("%d",r);//输出秘钥  
}
```

# 2013
## T1 计数问题-模拟，数位分离
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    int n, x, cnt = 0;  
    cin >> n >> x;  
    for (int i = 1; i <= n; i++)  
    {  
        int k = i;  
        while (k)  
        {  
            if (k % 10 == x) cnt++;  
            k /= 10;  
        }  
    }  
    cout << cnt;  
        return 0;  
}
```

## T2 表达式求值-模拟
本题可用标准双栈模板进行表达式求值，但是由于只涉及到加法和乘法且没有括号，所以可以直接进行模拟计算。
（1）直接模拟计算
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    long long ans = 0, t, num;  
    char ch;  
    cin >> t;  
    while(cin >> ch >> num){  
        if(ch == '+'){  
            ans = (ans+t) % 10000;  
            t = num;  
        }  
        else t = (t*num) % 10000;  
    }  
    ans = (ans+t) % 10000;  
    cout << ans;  
    return 0;  
}
```
（2）利用栈求中缀表达式的值
```cpp
/*author:闫晟淏*/
#include<bits/stdc++.h>
using namespace std;
stack<int> num;
stack<char> op;
void eval() {
    int b = num.top(); num.pop();
    int a = num.top(); num.pop();
    char c = op.top(); op.pop();
    int x;
    if (c == '+') x = (a + b)%10000;
    else if (c == '-') x = a - b;
    else if (c == '*') x = (a * b)%10000;
    else x = a / b;
    num.push(x);
}
unordered_map<char, int> pr{{'+',1},{'-',1},{'*',2},{'/',2}};
int main() {
    string str;
    cin >> str;
    for (int i = 0; i < str.size(); i++) {
        auto c = str[i];
        if (isdigit(c)) {
            int x = 0, j = i;
            while (j < str.size() && isdigit(str[j]))
            {
                x = x * 10 + str[j] - '0';
                j++;
            }
            i = j - 1;
            x=x%10000;
            num.push(x);
        }
        else if (c == '(') op.push(c);
        else if (c == ')') {
            while (op.top() != '(') eval();
            op.pop();
        }
        else {
            while (op.size() && op.top() != '(' && pr[c] <= pr[op.top()])
                eval();
            op.push(c);
        }
    }
    while (op.size()) eval();
    cout << num.top()%10000;
    return 0;
}
```

# 2014
## T1 珠心算测验-枚举，哈希表
枚举所有两两数之和，用哈希表（桶数组）标记，再枚举原数组检测是否在哈希表中。
（1）哈希表和动态数组写法
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	int n;
	cin >> n;
	vector<int> v(n+1);
	unordered_map<int, int> mp;
	for (int i = 1; i <= n; i++) cin >> v[i];
	for (int i = 1; i <= n; i++)
		for (int j = i+1; j <= n; j++)
			mp[v[i]+v[j]] = 1;
			
	int ans = 0;
	for (int i = 1; i <= n; i++)
	{
		if (mp[v[i]]) ans ++;
	}
	cout << ans;
	return 0;
}
```
（2）桶数组和静态数组写法
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int a[110], st[20010];  
  
int main()  
{  
    int n;  
    cin >> n;  
    for (int i = 1; i <= n; i++) cin >> a[i];  
    for (int i = 1; i <= n; i++)  
        for (int j = i+1; j <= n; j++)  
            st[a[i]+a[j]] = 1;  
    int ans = 0;  
    for (int i = 1; i <= n; i++)  
    {  
        if (st[a[i]]) ans ++;  
    }  
    cout << ans;  
    return 0;  
}
```
## T2 比例简化-枚举
L的范围在100内，枚举A'和B'的值，按题面要求模拟即可。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int gcd(int a, int b)  
{  
    if (a % b == 0) return b;  
    return gcd(b, a % b);  
}  
  
int main(){  
    int A, B, L;  
    cin >> A >> B >> L;  
    double minx = 1e9;  
    int a, b;  
    for (int i = 1; i <= L; i++)  
    {  
        for (int j = 1; j <= L; j++)  
        {  
            if (gcd(i, j) == 1)  
            {  
                double adivb = i*1.0/j;  
                double AdivB = A*1.0/B;  
                if (adivb >= AdivB && adivb - AdivB < minx)  
                    a = i, b = j, minx = adivb - AdivB;  
            }  
        }  
    }  
    cout << a << " " << b;  
    return 0;  
}
```

# 2015
## T1 金币-模拟
模拟金币发放规则求和
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    freopen("coin.in", "r",stdin);  
    freopen("coin.out", "w", stdout);  
    int n;  
    cin >> n;  
    int sum = 0, day = 0, q = 1; // 总，天数，钱数  
    for (int i = 1; i <= n; i++)  
    {  
        sum += q;  
        day++;  
        if (day == q) day = 0, q++;  
    }  
    cout << sum;  
        return 0;  
}
```

## T2 扫雷游戏-枚举，模拟
依次枚举每个格子，统计8个方向地雷数
```cpp
#include<bits/stdc++.h>  
using namespace std;  
const int N = 110;  
  
int n, m;  
char g[N][N];  
  
int main()  
{  
    scanf("%d%d", &n, &m);  
    for (int i = 0; i < n; i ++ ) cin >> g[i];  
        for (int i = 0; i < n; i ++ )  
    {  
        for (int j = 0; j < m; j ++ )  
            if (g[i][j] == '*') cout << '*';  
        else  
        {  
            int s = 0;  
            // 8方向枚举统计地雷数  
            for (int x = i - 1; x <= i + 1; x ++ )  
                for (int y = j - 1; y <= j + 1; y ++ )  
                    if (x != i || y != j)  
                    {  
                        if (x>=0&&x<n&&y>=0&&y<m&&g[x][y]=='*')                             s ++ ;  
                    }  
            cout << s;  
        }  
        cout << endl;  
    }  
        return 0;  
}
```

# 2016
## T1 买铅笔-模拟
由于只买同一种包装的铅笔，只需要对比所有铅笔的购买方案取最低花费即可
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    int n;  
    cin >> n;  
    int ans = 1e9;  
    for (int i = 1; i <= 3; i++)  
    {  
        int a, b; cin >> a >> b;  
        int cost = (n+a-1)/a*b;  // (a+b-1)/b可求a/b上取整结果  
        ans = min(ans, cost);  
    }  
    cout << ans;  
    return 0;  
}
```
## T2 回文日期-枚举，构造
只有8位数，且回文，可以枚举左边4位，右边4位构造出来，这样只需要枚举0~9999共一万个数，然后判断：
- 8位数日期是否合法
- 是否在范围内
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int months[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};  
  
bool check(int date) // 检测date是否合法  
{  
    int year = date / 10000;  
    int month = date % 10000 / 100;  
    int day = date % 100;  
    if (!month || month >= 13 || !day) return false;  
    if (month != 2 && day > months[month]) return false;  
    if (month == 2) // 2月特判  
    {  
        bool leap = year % 4 == 0 && year % 100 || year % 400 == 0;  
        if (day > 28 + leap) return false;  
    }  
    return true;  
}  
  
int main()  
{  
    int date1, date2;  
    cin >> date1 >> date2;  
    int res = 0;  
    for (int i = 0; i < 10000; i ++ )  
    {  
        int x = i, r = i;  
        for (int j = 0; j < 4; j ++ ) // 构造回文日期r  
        {  
            r = r * 10 + x % 10;  
            x /= 10;  
        }  
        if (r >= date1 && r <= date2 && check(r))             res ++ ;  
    }  
    cout << res;  
    return 0;  
}
```

# 2017
## T1 成绩-模拟
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    int a, b, c;  
    cin >> a >> b >> c;  
    cout << a*0.2+b*0.3+c*0.5;  
    return 0;  
}
```

## T2 图书管理员
对于每个需求码，查找包含它的最小图书编码，转换字符串处理，通过substr进行子串查找。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
const int N = 1010, INF = 1e8;  
int n, m;  
int id[N];  
  
int main()  
{  
    cin >> n >> m;  
    for (int i = 0; i < n; i ++ ) cin >> id[i];  
        while (m -- )  
    {  
        int l;  
        string num;  
        cin >> l >> num;  
        int res = INF;  
        for (int i = 0; i < n; i ++ )  
        {  
            string s = to_string(id[i]);  
            if (s.size() >= l && s.substr(s.size() - l) == num)  
                res = min(res, id[i]);  
        }  
        if (res == INF) res = -1;  
        cout << res << endl;  
    }  
    return 0;  
}
```
# 2018
## T1 标题统计-模拟，字符串
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	string s;
	getline(cin, s);
	int cnt = 0;
	for (int i = 0; i < s.size(); i++)
		if (s[i] != ' ' && s[i] != '\n')
			cnt ++;
	cout << cnt;
	return 0;
}
```

## T2 龙虎斗-模拟，枚举
经典模拟题，按照题面意思进行求解即可，详情见注释
```cpp
#include <bits/stdc++.h>  
#define int long long  
using namespace std;  
  
const int N=1e5+10;  
int a[N];  
  
signed main(){  
    freopen("fight.in","r",stdin);  
    freopen("fight.out","w",stdout);  
    int n;  
    cin >>n;  
    for(int i=1;i<=n;i++){  
        cin >>a[i];  
    }  
    int m,p1,s1,s2;  
    cin >>m>>p1>>s1>>s2;  
    a[p1]+=s1; // p1加入s1兵力  
    int l=0,r=0;  
    // 左边气势l  
    for(int i=1;i<m;i++){  
        l+=(m-i)*a[i];  
    }  
    // 右边气势r  
    for(int i=m+1;i<=n;i++){  
        r+=(i-m)*a[i];  
    }  
    int min_n=2e9;  
    int pos=0;  
    // 枚举所有位置，求加入s2后的最小气势差的position  
    for(int i=1;i<=n;i++){  
        int res_l=l,res_r=r;  
        if(i<m){  
            res_l+=(m-i)*s2;  
        }else{  
            res_r+=(i-m)*s2;  
        }  
        if(abs(res_l-res_r)<min_n){  
            pos=i;  
            min_n=abs(res_l-res_r);  
        }  
    }  
    cout <<pos;  
    return 0;  
}
```

## T3 摆渡车


# 2019
## T1 数字游戏-字符串
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    string s;  
    cin >> s;  
    cout << count(s.begin(), s.end(), '0');  
    return 0;  
}
```
## T2 公交换乘

从前往后扫描每条记录，同时用一个队列维护当前车次可以使用的优惠券区间：
- 如果当前记录是火车，则加入维护的优惠券区间；
- 如果当前记录是公交车，则线性扫描一遍队列中所有优惠券，找到第一个未被使用过的且大于等于当前价格的优惠券即可；
可以用一个bool数组对优惠券判重，以保证每张优惠券最多只被用一次。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
const int N = 100010;  
int n;  
struct Ticket  
{  
    int time, price; // 时间，票价  
}q[N];  
bool st[N];  
  
int main()  
{  
    scanf("%d", &n);  
        int res = 0;  
    for (int i = 0, l = 0, r = 0; i < n; i ++ )  
    {  
        int type, price, time; // 类型，价格，时间  
        scanf("%d%d%d", &type, &price, &time);  
        if (type == 0)  
        {  
            q[r ++ ] = {time, price};  
            res += price;  
        }  
        else  
        {  
            while (l < r && time - q[l].time > 45) l ++ ;  
                        bool success = false;  
            for (int j = l; j < r; j ++ )  
                if (!st[j] && q[j].price >= price)  
                {  
                    st[j] = true;  
                    success = true;  
                    break;  
                }  
            if (!success) res += price;  
        }  
    }  
        printf("%d\n", res);  
        return 0;  
}
```

## T3 纪念品
每天都可以无限买卖，只要保证今天买入的纪念品在明天卖出不亏本且获得最大利润。

初始金币：M
目标：让后面每一天的金币都尽可能变多

设第i天拥有的金币为$M_i$，怎么求第i+1天的金币$M_{i+1}$呢？
物品可以无限买，如果进去足够，且这个物品今天买明天卖能获利，那肯定一直买。

问题转换：在总钱数为$M_i$的情况下，对N个纪念品今天买，明天卖，最多获利多少钱。

**完全背包模型**
容量：$M_i$
体积：$price[i]$
价值：$price[i+1] - price[i]$

$M_{i+1} = M_i+$最大收益

时间复杂度：$O(TNM)$

```cpp
#include<bits/stdc++.h>
using namespace std;

constexpr int N = 110, M = 10010;

int t, n, m;
int price[N][N];
int dp[M];

int main(){
	cin >> t >> n >> m;
	for (int i = 1; i <= t; i ++)
		for (int j = 1; j <= n; j ++)
			cin >> price[i][j]; // 第i天的第j个物品价格
	
	for (int i = 1; i < t; i ++){
		// 因为我们是从第i天推到第i + 1天 所以i只能循环到t - 1
		memset(dp, 0, sizeof(dp));
		for (int j = 1; j <= n; j ++){ // 循环物品
			if (price[i + 1][j] > price[i][j]) // 优化
				for (int k = price[i][j]; k <= m; k ++)
					dp[k] = max(dp[k], dp[k - price[i][j]] + price[i + 1][j] - price[i][j]);
		} 
		m += dp[m];
	}
	cout << m << endl;
	return 0;
}
```
## T4 加工零件
题意：
给定一张包含n个点，m条边的无向图，再给定q个询问：ai, Li，判断是否存在一条从1号点走到ai号点恰好经过Li条边的路径。

关键思考：如果存在一条长度是L的路径，L > 0，那我们可以在任意一条边上来回走，就可以构造出L+2，L+4，L+6，...的路径。

因此当我们向判断是否存在长度是L的路径时，只需要判断是否存在长度小于L且奇偶性和L相同的路径即可。

因此问题转换成求：1号点到任意点长度为奇数的最短路和长度为偶数的最短路。

如何来分别求奇偶最短路呢？设1->v的路径长度是奇数，则再连一条边奇偶性会发生翻转。我们用`dist[k][0]`表示从1号点到k号点且长度为偶数的最短路径，`dist[k][1]`则表示长度为奇数的最短路径，每走到下一个点，奇偶性翻转即可。

```cpp
#include <bits/stdc++.h>
using namespace std;
typedef pair<int, int> PII;
const int N = 100010;
int n, m, query;
vector<int> adj[N];
int dist[N][2];

void bfs()
{
	memset(dist, 0x3f, sizeof dist); // 初始化dist为极大值
	queue<PII> q;
	q.push({1, 0}); // 起点1，阶段奇偶性0 
	dist[1][0] = 0; // 1->1，距离0
	while (!q.empty())
	{
		PII t = q.front();  
		q.pop();
		int ver = t.first; // 当前节点编号
		int type = t.second; // 奇偶类型
		
		for (int i = 0; i < adj[ver].size(); i++) // 遍历所有邻点
		{
			int j = adj[ver][i]; // 邻点 j
			if (dist[j][type ^ 1] > dist[ver][type] + 1)
			{
				dist[j][type ^ 1] = dist[ver][type] + 1;
				q.push({j, type ^ 1});
			}
		}
	}
}

int main()
{
	cin >> n >> m >> query;
	for (int i = 0; i < m; i++)
	{
		int a, b;
		cin >> a >> b;
		adj[a].push_back(b); // 无向图
		adj[b].push_back(a);
	}
	bfs();
	while (query--)
	{
		int a, L;
		cin >> a >> L;
		if (a == 1 && adj[1].empty()) // a为1且1没有邻居，不需要提供材料
			cout << "No\n";
		else if (dist[a][L & 1] <= L) // L奇数&1->1，L偶数&1->0
			cout << "Yes\n";
		else
			cout << "No\n";
	}
	return 0;
}

```

# 2020
## T1 优秀的拆分-二进制，位运算
只要是偶数就能进行拆分，因为偶数的二进制表示为`xxxxxxxx0`，`x`可为1或0.
所以有如下思路：n的范围是$10^7$，我们直接枚举0~30进行右移操作，n>>k&1 可以判断n的第k位是否为1，是的话输出1<<k表示这一位在10进制下的贡献。
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	int n;
	cin >> n;
	if (n & 1)
	{
		cout << -1;
		return 0;
	}
	for (int i = 30; i >= 0; i--)
		if (n >> i & 1)
			cout << (1 << i) << " ";
	return 0;
}
```
## T2 直播获奖-排序，桶
考虑到本题的数据范围，直接模拟的话每加入一个成绩都要重新维护序列的有序性，显然不行，仔细一看成绩的值最多600，此时应该想到桶数组进行有序性的维护。用桶数组记录每个成绩的人数，计算出获奖人数，从高到低去桶数组中查询分数即可。
```cpp
#include <bits/stdc++.h>
#define int long long
using namespace std;
int mp[610];
int get(int k){
    if(k==0) k=1;
    for(int i=600;i>=0;i--){
        k-=mp[i];
        if(k<=0) return i;
    }
    return -1;
}
signed main(){
    int n,w;
    cin >> n >> w;
    for(int i=1;i<=n;i++){
        int x;
        cin >> x;
        mp[x]++;
        cout << get(i*w/100) << " ";
    }
}
```
## T3 表达式
后缀表达式建树：
```cpp
#include <bits/stdc++.h>  
using namespace std;  
  
const int N = 1e5 + 10;  
vector<int> tree[N];  
char opt[N];     // 运算符  
int val[N];      // 叶子节点的数值  
bool isLeaf[N];  // 是否是叶子节点  
int idx = 0; // 节点编号计数  
stack<int> stk;  
  
// DFS 计算表达式值  
int dfs(int u) {  
    if (isLeaf[u]) return val[u]; // 叶子直接返回数值  
    int L = tree[u][0];  
    int R = tree[u][1];  
    if (opt[u] == '+') return dfs(L) + dfs(R);  
    if (opt[u] == '-') return dfs(L) - dfs(R);  
    if (opt[u] == '*') return dfs(L) * dfs(R);  
    if (opt[u] == '/') return dfs(L) / dfs(R);  
}  
  
int main() {  
    string s = "3 4 + 2 * 6 +"; // 后缀表达式  
    for (int i = 0; i < s.size(); i++)  
    {  
        if (isdigit(s[i])) // 数字  
        {             int sum = 0;  
            while (isdigit(s[i]))  
            {  
                sum = sum * 10 + s[i] - '0';  
                i++;  
            }  
            int id = ++idx;  
            val[id] = sum;  
            isLeaf[id] = true;  
            stk.push(id);  
        }         else if (s[i] != ' ') // 运算符  
        {  
            int b = stk.top(); stk.pop();  
            int a = stk.top(); stk.pop();  
            int id = ++idx;  
            opt[id] = s[i];  
            isLeaf[id] = false;  
            tree[id].push_back(a);  
            tree[id].push_back(b);  
            stk.push(id);  
        }  
    }  
        int root = stk.top();  
    cout << "表达式值 = " << dfs(root) << "\n";  
    return 0;  
}
```

【参考代码】
```cpp
#include <bits/stdc++.h>
using namespace std;

const int N = 1e6 + 10;
vector<int> tree[N];
char opt[N];     // 运算符
int val[N], w[N], vis[N];      // 叶子节点的数值
bool isLeaf[N];  // 是否是叶子节点
int idx = 0; 	// 节点编号计数
stack<int> stk;

// dfs 计算表达式值
int dfs(int u) {
	if (isLeaf[u]) return val[u]; // 叶子直接返回数值
	if (opt[u] == '!')
	{
		int e = tree[u][0];
		return w[u] = !dfs(e);
	}
	else if (opt[u] == '&')
	{
		int t = 1;
		int a = tree[u][0], b = tree[u][1];
		t = t & dfs(a) & dfs(b);
		return w[u] = t;
	}
	else if (opt[u] == '|')
	{
		int t = 0;
		int a = tree[u][0], b = tree[u][1];
		t = t | dfs(a) | dfs(b);
		return w[u] = t;
	}
}

void dfs2(int u)
{
	if (isLeaf[u])
	{
		vis[u] = 1;
		return;
	}
	if (opt[u] == '!')
	{
		int e = tree[u][0];
		dfs2(e);
	}
	else if (opt[u] == '&')
	{
		int a = tree[u][0];
		int b = tree[u][1];
		if (w[a] == 1) dfs2(b);
		if (w[b] == 1) dfs2(a);
	}
	else if (opt[u] == '|')
	{
		int a = tree[u][0];
		int b = tree[u][1];
		if (w[a] == 0) dfs2(b);
		if (w[b] == 0) dfs2(a);
	}
}

int main() {
	string s;
	getline(cin, s);
	int n; cin >> n;
	idx = n;
	for (int i = 1; i <= n; i++) cin >> val[i], w[i] = val[i];
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] == ' ') continue;
		if (s[i] == 'x') // 数字
		{ 
			i++;
			int id = 0;
			while (isdigit(s[i]))
			{
				id = id * 10 + s[i] - '0';
				i++;
			}
			isLeaf[id] = true;
			stk.push(id); // val[id]取到值
		} 
		else if (s[i]  == '!') // 运算符
		{
			int a = stk.top(); stk.pop();
			opt[++idx] = s[i];
			isLeaf[idx] = false;
			tree[idx].push_back(a);
			stk.push(idx);
		}
		else if (s[i] == '|' || s[i] == '&')
		{
			int b = stk.top(); stk.pop();
			int a = stk.top(); stk.pop();
			opt[++idx] = s[i];
			tree[idx].push_back(a);
			tree[idx].push_back(b);
			stk.push(idx);
		}
	}
	
	int root = stk.top();
	int ans = dfs(root);
	dfs2(root);
	int q;
	cin >> q;
	while (q--)
	{
		int x; cin >> x;
		if (vis[x]) cout << !ans << endl;
		else cout << ans << endl;
	}
	return 0;
}

```


## T4 方格取数
思路：

首先考虑最普通的方格取数问题进行DP，转移方程有：

`f[i][j] = max(f[i][j-1], f[i-1][j], f[i+1][j])`

但是这样做是有后效性的，因为更新完`f[i+1][j]`后可能会重复更新`f[i][j]`。

鉴于此，我们要考虑新的做法，即把状态进一步拆分：

令：

- `f[i][j][0]` 表示走到`(i,j)`，且从左边走过来的路径最大值
- `f[i][j][1]` 表示走到`(i,j)`，且从上面走过来的路径最大值
- `f[i][j][2]` 表示走到`(i,j)`，且从下面走过来的路径最大值

对于`f[i][j][0]`，显然是由` (i,j-1) `这个点向右走了一格而来，所以有：

`f[i][j][0] = max(f[i][j-1][0], f[i][j-1][1], f[i][j-1][2])`

对于` f[i][j][1]`，显然是由 `(i-1,j)` 这个点向下走了一格而来，所以有：

`f[i][j][1] = max(f[i-1][j][0], f[i-1][j][1])`，注意并不能从`f[i-1][j][2]`转移！

对于` f[i][j][2]`，显然是由` (i+1,j)` 这个点向上走了一格而来，所以有：

`f[i][j][2] = max(f[i+1][j][0], f[i+1][j][2])`，注意并不能从`f[i+1][j][1]`转移！

时间复杂度：`O(3nm)`

基础情况：

由于是求最大值，可以初始化`f[N][N][3]`整个数组为极小值，注意开 long long。
`f[1][1][0] = g[1][1]`
`f[1][1][1] = g[1][1]`
```cpp
#include <bits/stdc++.h>
#define ll long long
#define LINF 0x3f3f3f3f3f3f3f3f
using namespace std;
const int N = 1e3+10 ;
const ll MIN = -LINF;
int g[N][N];
ll f[N][N][3];//由  右 上 下 转移而来 注意开long long
int n,m;

int main() 
{
	cin >> n >> m;	
	for(int i=1;i<=n;i++)
		for(int j=1;j<=m;j++)
			cin >> g[i][j];
	//初始化
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			for(int k=0;k<3;k++)
				f[i][j][k] = MIN;
	
	f[1][1][0] = (ll)g[1][1];
	f[1][1][1] = (ll)g[1][1];
	f[1][1][2] = (ll)g[1][1];
	
	for (int i = 2; i <= n; i++)
		f[i][1][2] = f[i-1][1][2] + g[i][1];
	
	for (int j = 2; j <= m; j ++ )
	{
		for (int i = 1; i <= n; i ++)
		{
			f[i][j][0] = max({f[i][j-1][0],f[i][j-1][1],f[i][j-1][2]}) + (ll)g[i][j];
			f[i][j][1] = max(f[i-1][j][0],f[i-1][j][1]) + (ll)g[i][j];
		}
		for (int i = n-1; i >= 0; i --)
		{
			f[i][j][2] = max(f[i+1][j][0],f[i+1][j][2]) + (ll)g[i][j];
		}
	}
	cout << max(f[n][m][0], f[n][m][1]);
	return 0;
}

```

# 2021
## T1 分糖果-数学
本题是一道数学思维题，认真画图思考后可得到如下性质：
- 若L/n == R/n 则说明L~R之间最多n个，此时可获得最多奖励为：R%n
- 否则一定可以获得 n-1个
```cpp
#include<bits/stdc++.h>
using namespace std;

int main()
{
	freopen("candy.in", "r", stdin);
	freopen("candy.out", "w", stdout);
	int n, L, R;
	cin >> n >> L >> R;
	if (L / n == R / n) // l~r不够n个
		cout << R % n;
	else
		cout << n - 1;
	
	return 0;
}
```

## T2 插入排序-模拟，排序
本题的难点在于卡时间复杂度，需要我们用类似冒泡排序的方法对新加入的数在线性时间复杂度内维护有序性，涉及到多关键字的排序，除了交换值之外，还需要把对应的下标映射一并进行交换，这样才能在修改操作中正确找到目标。
```cpp
#include<bits/stdc++.h>
#define io(x); freopen(x".in", "r", stdin), freopen(x".out", "w", stdout);

using namespace std;

const int N = 8010;

int n, m;
int a[N], b[N], c[N];

bool cmp(int x, int y)
{
	if (a[x] != a[y]) return a[x] < a[y];
	return x < y;
}

void upd(int k)
{
	while (k > 1 && (a[b[k - 1]] > a[b[k]] || a[b[k - 1]] == a[b[k]] && b[k - 1] > b[k]))
	{
		swap(b[k - 1], b[k]);
		swap(c[b[k - 1]], c[b[k]]);
		k -- ;
	}
	while (k < n && (a[b[k + 1]] < a[b[k]] || a[b[k + 1]] == a[b[k]] && b[k + 1] < b[k]))
	{
		swap(b[k + 1], b[k]);
		swap(c[b[k + 1]], c[b[k]]);
		k ++ ;
	}
}

int main()
{
	io("sort");
	cin >> n >> m;
	for (int i = 1; i <= n; i ++ )
	{
		cin >> a[i];
		b[i] = i; // 存a[i]的下标
	}
	
	sort(b + 1, b + n + 1, cmp); 
	for (int i = 1; i <= n; i ++ )
		c[b[i]] = i; // 存原数组第 k 个位置的数在排完序后的位置
	
	while (m -- )
	{
		int op, x, v;
		cin >> op >> x;
		if (op == 1)
		{
			cin >> v;
			a[x] = v; // 修改
			int k = c[x]; // 获取排完序后的位置
			upd(k);
		}
		else cout << c[x] << endl;
	}
	
	return 0;
}
```

## T3 网络连接
所有的非法情况：
前导零，如19.19.08.10:0或11.4.5.14:00;
字符位置不对，如19.19:8.10.0
多余字符，如19.19.8.10:0:
数字大小超出范围，如19.19.8.10:114514或1919.8.1.0:11
```cpp
#include<bits/stdc++.h>  
using namespace std;  
  
bool check(string addr)  
{  
    vector<int> w(5, -1);  
    sscanf(addr.c_str(), "%d.%d.%d.%d:%d", &w[0], &w[1], &w[2], &w[3], &w[4]);  
    for (int i = 0; i < 5; i ++ )  
    {  
        if (w[i] < 0) return false;  
        if (i < 4 && w[i] > 255) return false;  
        if (i == 4 && w[i] > 65535) return false;  
    }  
    char str[100];  
    sprintf(str, "%d.%d.%d.%d:%d", w[0], w[1], w[2], w[3], w[4]);  
    return str == addr;  
}  
  
int main()  
{  
    int n;  
    cin >> n;  
    unordered_map<string, int> hash;  
    for (int i = 1; i <= n; i ++ )  
    {  
        string type, addr;  
        cin >> type >> addr;  
        if (type == "Server")  
        {  
            if (!check(addr)) puts("ERR");  
            else if (hash.count(addr)) puts("FAIL");  
            else  
            {  
                hash[addr] = i;  
                puts("OK");  
            }  
        }  
        else  
        {  
            if (!check(addr)) puts("ERR");  
            else if (!hash.count(addr)) puts("FAIL");  
            else cout << hash[addr] << endl;  
        }  
    }  
        return 0;  
}
```
## T4 小熊的果篮
```cpp
#include<bits/stdc++.h>  
using namespace std;  
const int N = 2e5 + 10;  
int n;  
int w[N];  
int l[N],r[N],bl[N],br[N];  
  
int main()  
{  
    freopen("bear.in","r",stdin);  
    freopen("bear.out","w",stdout);  
    cin >> n;  
    for (int i = 1; i <= n; i++) cin >> w[i];  
    int cnt = 0; // 块数量  
    for (int i = 1, last = n; i <= n; i++)  
    {  
        int j = i + 1;  
        while (j <= n && w[j] == w[i]) j++; // i ~ j-1相同  
        // 建立i~j-1的双向链表  
        for (int k = i; k <= j-1; k++)  
        {  
            l[k] = k-1;  
            r[k] = k+1;  
        }  
        l[i] = j-1, r[j-1] = i; // 首位相连，双向循环链表  
        // 块链表，每个块的最后一个位置作为标志，last上一个块，j-1当前块  
        br[last] = j-1, bl[j-1] = last;  
        i = j-1; // 跳到下一个块  
        last = j-1; // 更新上一个块的指针  
        cnt ++;  
    }  
    int head = br[n]; // 最后一个块的右边 = 第一个块的末尾指针  
    while (true)  
    {  
        int new_head = -1;  
        // 遍历所有块进行输出与合并  
        for (int t = cnt, i = head; t != 0; t--, i = br[i])  
        {  
            int k = r[i]; // 块内的第一个元素  
            cout << k << " ";  
            if (k != i) // 第一个不等于最后一个，说明块内元素个数>1  
            {  
                // 删除输出了的元素  
                r[l[k]] = r[k], l[r[k]] = l[k];  
                if (new_head == -1) new_head = i;  
            }  
            else // 块内只有一个元素，整块删除  
            {  
                br[bl[k]] = br[k], bl[br[k]] = bl[k];  
                cnt --;  
            }  
        }  
        if (cnt == 0) break;  
                // 合并块  
        head = new_head;  
        if (cnt > 1)  
        {  
            new_head = -1;  
            for (int t = cnt, i = head; t > 1; t--, i = br[i])  
            {  
                int j = br[i];  
                if (w[i]==w[j])  
                {  
                    int ri = r[i], rj = r[j];                     // i合并到j  
                    r[i] = rj, l[rj] = i;  
                    l[ri] = j, r[j] = ri;  
                    // 删除i  
                    br[bl[i]] = br[i];  
                    bl[br[i]] = bl[i];  
                    cnt--;  
                }  
                else if (new_head == -1) new_head = i;  
            }  
            if (new_head == -1) new_head = bl[head];  
            head = new_head;  
        }  
        cout << endl;  
    }  
        return 0;  
}
```

# 2022
## T1 乘方-模拟，数学
注意特殊情况：a = 1，其他正常模拟即可
```cpp
#include <bits/stdc++.h>  
using namespace std;  
  
int main()  
{  
    int a,b;  
    cin >>a>>b;  
    if(a==1)  
    {  
        cout <<"1";  
        return 0;  
    }  
    int ans=1;  
    for(int i=1;i<=b;i++)  
    {  
        ans*=a;  
        if(ans>1e9)  
        {  
            cout <<"-1";  
            return 0;  
        }  
    }  
    cout <<ans;  
    return 0;  
}
```

## T2 解密
推导过程如下：

$ed = (p-1)(q-1)+1$
$ed = pq-p-q+2$

$∵ n = pq$
$∴ ed = n-p-q+2$

已知n、e、d，可知p+q为：
$ed = n-p-q+2$
$p+q+ed = n+2$
$p+q = n-ed+2$

因此：
$pq = n$
$p+q = n-ed+2$

用**韦达定理**可以求得p和q的值，也可以直接推导求解：

令$n-ed+2 = m$
$pq = n$
$p+q = m$

完全平方公式：
> $(a+b)^2 = a^2 + 2ab + b^2$
> $(a-b)^2 = a^2 - 2ab - b^2$

所以有：
$(p+q)^2=p^2-2pq+q^2+4pq$ 得：
$p^2-2pq+q^2 = m^2 - 4pq$
$(p-q)^2 = m^2 - 4pq$
$p-q = sqrt(m^2 - 4pq)$

令 $sqrt(m^2-4pq)  = t$ 得：
$(1)p+q = m$
$(2)p-q = t$
解得：
$p = (m+t)/2$
$q = (m-t)/2$

```cpp
#include<bits/stdc++.h>  
using namespace std;  
typedef long long ll;  
int k;  
ll n, d, e;  
int main()  
{  
    cin >> k;  
    while(k -- ){  
        cin >> n >> d >> e;  
        ll m = n - e * d + 2;  
        ll t = sqrt(m * m - 4 * n);  
        ll p = (m - t) / 2;  
        ll q = (t + m) / 2;  
        if (p > 0 && q > 0 && p + q == m)  
            cout << p << ' ' << q << endl;  
        else  
            cout << "NO" << endl;  
    }  
}
```

## T3 逻辑表达式
本题统计的是表达式中&和|两种运算符被短路的次数，假设分别为x和y，对于表达式中的某个非叶子节点（一定是两种运算符中的某种），假设其左右儿子分别为𝑙和𝑟，分类如下：
根是&：
左儿子𝑙是0：此时以该节点为根的子树构成的&子表达式被短路，右子树被忽略，根的短路次数为左儿子短路次数加1：{l.x + 1, l.y}；
左儿子𝑙是1：此时子表达式未短路，根的短路次数为左右儿子短路次数之和：{l.x + r.x, l.y + r.y}；

根是 |：
左儿子𝑙是1：此时|子表达式被短路，右子树被忽略，根的短路次数为左儿子短路次数加1：{l.x, l,y + 1}；
左儿子𝑙是0：此时子表达式未短路，根的短路次数为左儿右儿子短路次数之和：{l.x + r.x, l.y + r.y}；

```cpp
#include<bits/stdc++.h>
#define x first     
#define y second

using namespace std;
typedef pair<int, int> PII;

stack<int> lg;              // 值栈（logic value stack）：存储子表达式的布尔值（0 或 1）
stack<char> op;              // 运算符栈：存储 '('、'&'、'|'
stack<PII> cnt;              // 计数栈：与值栈同步，存储 {&短路次数, |短路次数}
unordered_map<char, int> pri{{'(', 0}, {'|', 1}, {'&', 2}}; 
// 运算符优先级表：括号最低，| 次之，& 最高（保证 & 比 | 先算）

// 对栈顶的一个运算进行计算，并更新计数
void eval() {
    // 按照“值栈”的出栈顺序，先取右操作数 a，再取左操作数 b
    int a = lg.top(); lg.pop();       // 右操作数的值
    PII acnt = cnt.top(); cnt.pop();  // 右操作数的短路计数

    int b = lg.top(); lg.pop();       // 左操作数的值
    PII bcnt = cnt.top(); cnt.pop();  // 左操作数的短路计数

    char c = op.top(); op.pop();      // 运算符（'&' 或 '|'）

    if (c == '&') {
        // 计算 a & b 的布尔值（位运算，保证输入是 0/1）
        lg.push(a & b);

        // 短路规则：
        // 对于 b & a（左边是 b）：如果 b == 0，短路发生（右边不用算）
        // 因为 b 是左操作数，而栈是“右先出”，所以这里判断 b
        if (!b) {
            // 左为 0，短路一次（bcnt.x + 1），右边的计数不再累加
            cnt.push({bcnt.x + 1, bcnt.y});
        } else {
            // 左不为 0（即左为 1），要执行右边，所以两边计数相加
            cnt.push({acnt.x + bcnt.x, acnt.y + bcnt.y});
        }
    } else if (c == '|') {
        // 计算 a | b
        lg.push(a | b);

        // 短路规则：
        // 对于 b | a：如果 b == 1，短路发生（右边不用算）
        if (b == 1) {
            // 左为 1，短路一次（bcnt.y + 1），右边的计数不加
            cnt.push({bcnt.x, bcnt.y + 1});
        } else {
            // 左不为 1（即左为 0），要执行右边，所以两边计数相加
            cnt.push({acnt.x + bcnt.x, acnt.y + bcnt.y});
        }
    }
}

int main() {
    string expr;
    cin >> expr; // 输入一个布尔表达式，例如 "(1&0)|1"

    for (int i = 0; i < expr.size(); ++i) {
        if (isdigit(expr[i]))  { 
            // 如果是数字（0 或 1）
            // 压入值栈（转成 int）
            lg.push(expr[i] - '0');
            // 对应计数栈压入 {0,0}（叶子节点没有短路次数）
            cnt.push({0, 0});

        } else if (expr[i] == '(') {
            // 左括号直接压入运算符栈
            op.push(expr[i]);

        } else if (expr[i] == ')') {
            // 遇到右括号：不断计算，直到遇到左括号
            while (op.top() != '(') eval();
            op.pop(); // 弹出左括号

        } else { 
            // 运算符（'&' 或 '|'）
            // 当运算符栈顶的优先级 >= 当前运算符时，先执行栈顶的运算（保证左结合）
            while (op.size() && pri[op.top()] >= pri[expr[i]]) eval();
            // 再将当前运算符压入栈
            op.push(expr[i]);
        }
    }

    // 扫描完后，把剩余的运算符全部计算
    while (op.size()) eval();

    // 栈顶的值就是最终表达式的结果
    cout << lg.top() << '\n';
    // 栈顶的计数就是总的短路次数
    cout << cnt.top().x << " " << cnt.top().y << '\n';

    return 0;
}
```

## T4 上升点列-动态规划，LIS变形
```cpp
#include <bits/stdc++.h>
#define x first
#define y second
using ll = long long;
using namespace std;

const int N = 510;
int dp[N][N]; // dp[i][j]表示以i点结尾添加j个点的最大长度
pair<int,int> a[N];

int main() {
	freopen("point.in", "r", stdin);
	freopen("point.out", "w", stdout);
	int n, k;
	cin >> n >> k;
	for (int i = 1; i <= n; i++) 
		cin >> a[i].x >> a[i].y;
	sort(a+1, a+1+n);
	// 初始化
	for (int i = 1; i <= n; i++)
		for (int j = 0; j <= k; j++)
			dp[i][j] = j + 1; // 以i结尾添加j个点的长度至少j+1
	// 状态计算
	for (int i = 2; i <= n; i++)  // 从2号点开始枚举
		for (int j = 0; j < i; j++) // 枚举1~i-1号点
		{
			if (a[j].y > a[i].y) continue;
			int d = a[i].x - a[j].x + a[i].y - a[j].y - 1;
			for (int p = d; p <= k; p++) // 枚举添加的点数量
				dp[i][p] = max(dp[i][p], dp[j][p-d]+d+1);
		}
	// ans
	int ans = 0;
	for (int i = 1; i <= n; i++) ans = max(ans, dp[i][k]);
	cout << ans;
	return 0;
}
```
# 2023
## T1 小苹果-数学，模拟
通过模拟样例和计算可以发现，每次拿走的苹果数量为 n/3 上取整。循环模拟取苹果的过程，注意若`n%3==1`时会取走第n个苹果。
```cpp
#include<bits/stdc++.h>  
using namespace std;  
int n;  
int main(){  
    cin >> n;  
    int cnt = 0,ans = 0; // 分别对应答案1，答案2  
    while(n > 0){ // 直到所有苹果被取完  
        cnt ++; // 取苹果的轮数加1  
        int sum = ceil((double)n / 3); // 需要拿走的苹果数量  
        if(n % 3 == 1 && ans == 0) ans = cnt; // 如果此刻可以拿走编号为 n 的苹果，记录  
        n -= sum; // 取走苹果  
    }  
    cout << cnt << " " << ans;  
    return 0;  
}
```
## T2 公路-贪心，模拟
本题的贪心策略为：如果一个点买油的单价是最便宜的，那么我们走完后面全部路程的油都在这里买即可。也就是如果ai最小，则`[i,n]`的路我们都用这个油。此时`[1,i]`的路我们该怎么买油呢?也是一样的过程。按照如上流程模拟计算答案。
```cpp
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
int n,d;
const int N = 1e5 + 10;
ll ans=0;
ll v[N], a[N], y[N];
ll youjia;//油价
int main() 
{
	cin >> n >> d;
	for (int i=2;i<=n;i++){
		cin >> v[i];
		v[i]+=v[i-1];//计算站点i到起点的距离
		y[i]=ceil(1.0*v[i]/d);//从起点到当前站点的油耗
	}
	for (int i=1;i<=n;i++) cin >> a[i];
	
	youjia=a[1];//出发先在站点1买油
	for (int i=2;i<=n;i++){
		ans+=youjia*(y[i]-y[i-1]);//增加本段用油的油价
		youjia=min(youjia,a[i]);//如果当前站点油比用的油便宜，换成当前站点买
	}
	cout<<ans;
	return 0;
}
```
