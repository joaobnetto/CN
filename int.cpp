#include <bits/stdc++.h>
#include <boost/math/tools/polynomial.hpp>
#define INF INT_MAX
#define EPS 10e-9
#define IMG 1001

using namespace std;
using namespace boost::math;
using namespace boost::math::tools;
using boost::lexical_cast;

/*
*	Função auxiliar para imprimir o polinomio.
*/

void imprimePolinomio( const polynomial < long double > &res ){
	cout << "P(x)= ";
	for(int i = 0;i < res.size();++i){
		if(i == 0){
			cout << res[i];
			continue;
		}
		double x = res[i];
		bool sinal = x < 0;
		if(x < 0)
			x *= -1;
		if(i == 1)
			cout << (sinal ? " - " : " + ") << x << "x";
		else
			cout << (sinal ? " - " : " + ") << x << "x^" << i;
	}
	cout << endl;
}

/*
*	Função de Gauss para descobrir o resultado de um sistema linear.
*	Usando o do cp-algorithms.com
*/
int gauss (vector< vector<long double> > a, vector<long double> & ans) {
    int n = (int) a.size();
    int m = (int) a[0].size() - 1;

    vector<int> where (m, -1);
    for (int col=0, row=0; col<m && row<n; ++col) {
        int sel = row;
        for (int i=row; i<n; ++i)
            if (abs (a[i][col]) > abs (a[sel][col]))
                sel = i;
        if (abs (a[sel][col]) < EPS)
            continue;
        for (int i=col; i<=m; ++i)
            swap (a[sel][i], a[row][i]);
        where[col] = row;

        for (int i=0; i<n; ++i)
            if (i != row) {
                long double c = a[i][col] / a[row][col];
                for (int j=col; j<=m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        ++row;
    }

    ans.assign (m, 0);
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i]][m] / a[where[i]][i];
    for (int i=0; i<n; ++i) {
        long double sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * a[i][j];
        if (abs (sum - a[i][m]) > EPS)
            return 0;
    }

    for (int i=0; i<m; ++i)
        if (where[i] == -1)
            return INF;
    return 1;
}

/*
*	Função que monta um sistema linear baseado nos vetores dados
*	E resolve o sistema linear.
*	Retorna um vetor vazio se o sistema não existe.
*/
polynomial < long double > resolverGauss( const vector< long double > &a, const vector < long double > &b ){
	const int N = a.size();

	// Crio um novo sistema para ser resolvido
	vector< vector< long double > > sistema(N);
	vector< long double > res;

	for(int i = 0;i < N;++i){
		// x^0 = 1, logo tmp começa com 1
		long double tmp = 1;
		// Mudo o tamanho da linha pra N+1, pois tem a0x^0...anx^nn = bi
		sistema[i].resize(N+1);
		for(int j = 0;j < N;++j){
			// Coloco o x no lugar dele
			sistema[i][j] = tmp;
			// Multiplico o x atual por ele mesmo, pra conseguir as potencias
			tmp *= a[i];
		}
		// A ultima posição é f(x)
		sistema[i][N] = b[i];
	}
	// Chamo gauss pra resolver, e retorno dependendo do número de soluções.
	int gaussAns = gauss( sistema, res );
	if(gaussAns == 0){
		cout << "Não há soluções para esse sistema por Gauss" << endl; 
		res.clear();
	}
	polynomial< double > p(res.begin(), res.end());
	return p;
}

/*
*	Essa função calcula um Lk(x) de acordo com laGrange, ele tem como argumentos
*	o k, e o vetor com os X.
*	Estou reutilizando o vetor com os X normais obtidos na main.
*/

polynomial< long double > calcularLk( const int k, const vector< long double > &a ){
	/*
	*	Começo a parte de baixo com 1 e a parte de cima também com 1, pois 
	*	é o elemento neutro da multiplicação.
	*/
	long double bottom = 1;
	polynomial< long double > top = {1};
	for(int i = 0;i < a.size();++i){
		/* 
		*	Pela regra de LaGrange, multiplico todos menos quando i = k,
		*	pois se não, x - xk = 0;
		*/
		if(i == k)
			continue;
		// Multiplicando o de baixo, de acordo com a formula.
		bottom *= (a[k] - a[i]);
		/*
		*	Crio um polinomio em cima, x - a[i][1], 
		*	desse modo, posso multiplicar depois.
		*/
		polynomial< long double > tmp = {{-1*a[i], 1}};
		/*
		*	Multiplico o polinomio atual pelo que tenho até agora
		*/
		top *= tmp;
	}
	// No final é só dividir o polinomio pela parte de baixo.
	top /= bottom;
	return top;
}

/*
*	Essa função acha um polinomio de interpolação usando a fórmula de LaGrange
*	Tem como argumentos o vetor com os x, e o vetor com os f(x)
*/
polynomial< long double > resolverLaGrange( const vector< long double > &a, const vector < long double > &b ){
	// Começo com 0 pois é o elemento neutro da adição.
	polynomial< long double > p = {{0}};

	for(int i = 0;i < a.size();++i){
		// Calculo Li(x)
		polynomial< long double > tmp = calcularLk(i, a);
		// Multiplico ele pelo f(x).
		tmp *= b[i];
		// Adiciono ele ao vetor atual
		p += tmp;
	}
	return p;
}

/*
*	Função que calcula a tabela que é usada pelo metódo de Newton
*/
vector< vector< long double > >calcularTabela( const vector < long double > &a, const vector < long double > &b ){
	const int N = a.size();

	vector < vector < long double > > table(N);
	// Primeiramente colocamos a primeira coluna com os resultados de f(x)
	for(int i = 0;i < N;++i){
		table[i].push_back(b[i]);
	}

	// Depois calculamos usando a formula abaixo, pegando da coluna anterior
	// Linha atual e anterior.
	for(int i = 1;i < N;++i){
		int k = 0;
		for(int j = i;j < N;++j){
			table[j].push_back( (table[j][i-1] - table[j-1][i-1]) / (a[j] - a[k++]) );
		}
	}

	return table;
}

/* 
*	Metodo de Newton para resolver interpolação, constroi a tabela e acha o polinomio
*/
polynomial< long double > resolverNewton( const vector < long double > &a, const vector < long double > &b ){
	const int N = a.size();
	vector< vector< long double > > table = calcularTabela( a, b );
	/*
	*	Inicio o atual, que vai conter o polinomio até aqui (x-x0)(x-x1).. etc como 1, elemento neutro
	*	E res começa como 0, elemento neutro da adição, pois vai ser o polinomio resultado.
	*/
	polynomial< long double > atual = {{1}}, res = {{0}};
	for(int i = 0;i < N;++i){
		// Sempre calculo o atual vezes a tabela, crio um polinomio temporario e multiplico, obtendo o resultado
		res += (atual * table[i][i]);
		polynomial< long double > tmp = {{a[i]*-1, 1}};
		atual *= tmp;
	}
	return res;
}


int main(){
	int n;
	long double x, v;
	cout << "Digite o número de pontos:\n";
	cin >> n;
	vector < long double >a(n), b(n);
	cout << "Digite x e f(x) para cada ponto:\n";
	for(int i = 0;i < n;++i){
		cin >> x >> v;
		a[i] = x;
		b[i] = v;
	}

	polynomial< long double > poliGauss = resolverGauss( a, b );
	polynomial< long double > poliGrange = resolverLaGrange( a, b );
	polynomial< long double > poliNewton = resolverNewton( a, b );

	cout << "Metodo de Gauss\n";
	imprimePolinomio( poliGauss );

	cout << "Metodo de Grange\n";
	imprimePolinomio( poliGrange );

	cout << "Metodo de Newton\n";
	imprimePolinomio( poliNewton );
}