#include <vector>

#include <iomanip>

#include <sstream>

#include <iostream>

#include <string>

// Define DEBUG_MODE to enable or disable debug messages.
#define DEBUG_MODE 1

#if DEBUG_MODE
#define DEBUG_PRINT(x) std::cout << "     ->" << x << std::endl
#else
#define DEBUG_PRINT(x)
#endif

void print_spaces(int spaces) {
    for (int k = 0; k < spaces; k++)
        std::cout << " ";
}



//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */

class fraction {
    public: 
    int get_numerator() const;
    fraction(); // sets numerator to 0 and denominator to 1
    fraction(int); // sets numerator to the paramamter and denominator to 1
    fraction(int, int); //sets the numerator and denominator to the paramaters
    fraction(const std::string & input); // sets the fraction to a string
    void reduce_fraction(); // reduces fraction to its simplest form.The function divides the numerator and denomitaor by the GCD.
    int gcd_Euclidean(int, int) const; //helper function for the reduce_fraction. determinbes the Greatest Common Divisor of two integers.
    friend int length(const fraction &) ; //finds the length of the fraction using string stream; useful for printing the fraction in a matrix to calculate padding.
    friend std::string to_string(const fraction & );

    friend fraction operator + (const fraction & ,
        const fraction & );
    friend std::ostream & operator << (std::ostream & os,
        const fraction & f1);
    friend std::istream & operator >> (std::istream & is, fraction & f);
    friend fraction operator - (const fraction & ,
        const fraction & );
    friend fraction operator * (const fraction & ,
        const fraction & );
    friend fraction operator / (const fraction & ,
        const fraction & );
    friend fraction operator + (const fraction & , int);
    friend fraction operator - (const fraction & , int);
    friend fraction operator * (const fraction & , int);
    friend fraction operator / (const fraction & , int);
    friend void operator += (const fraction & ,
        const fraction & );
    friend void operator -= (const fraction & ,
        const fraction & );
    friend void operator *= (const fraction & ,
        const fraction & );
    friend void operator /= (const fraction & ,
        const fraction & );
    friend bool operator!(const fraction & );
    friend bool operator != (const fraction & ,
        const fraction & );
    friend bool operator != (const fraction & , int);
    friend bool operator == (const fraction & ,
        const fraction & );
    friend bool operator == (const fraction & , int);


    private: 
    int numerator;
    int denominator;
};

int fraction::get_numerator() const {
    return numerator;
}

fraction::fraction() {
    numerator = 0;
    denominator = 1;
}

fraction::fraction(int n): numerator(n), denominator(1) {}

fraction::fraction(int n, int d): numerator(n), denominator(d) {
    if (denominator == 0) {
        exit(4);
    }
    reduce_fraction();
}

//For a string. ex) 5/4
fraction::fraction(const std::string & input) {
    std::stringstream ss(input);
    char slash;
    if (ss >> numerator) {
        if (ss >> slash && slash == '/' && ss >> denominator) {
            if (denominator == 0) {
                exit(8);
            }
        } else {
            denominator = 1;
        }
    } else {
        exit(6);
    }
    reduce_fraction();
}

int fraction::gcd_Euclidean(int a, int b) const {
    while (b != 0) {
        int remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

void fraction::reduce_fraction() {
    int gcd_value = gcd_Euclidean(numerator, denominator);
    if (gcd_value != 1) {
        numerator /= gcd_value;
        denominator /= gcd_value;
    }
    if (!numerator) {
        denominator = 1;
    }
    if (denominator < 0) {
        //ensure only the numerator can be negative. The negative is transfered from the denominator to the numerator.
        numerator *= -1;
        denominator *= -1;
    }
}

int length(const fraction & f) {
    std::ostringstream oss;
    oss << f;  // Use the overloaded << operator
    return oss.str().length();  // Get the length of the resulting string
}


std::string to_string(const fraction & f) {
    return (std::to_string(f.numerator) + "/" + std::to_string(f.denominator));
}

std::istream & operator >> (std::istream & is, fraction & f) {
    std::string input;
    is >> input;
    std::stringstream ss(input);
    char slash;
    if (ss >> f.numerator) {
        if (ss >> slash && slash == '/' && ss >> f.denominator) {
            if (f.denominator == 0) {
                exit(0);
            }
        } else {
            f.denominator = 1;
        }
    } else {
        exit(6);
    }
    f.reduce_fraction();
    return is;
}

fraction operator + (const fraction & f1,
    const fraction & f2) {
    fraction f3;
    if (f1.denominator == f2.denominator) {
        f3.numerator = f1.numerator + f2.numerator;
        f3.denominator = f1.denominator;
    } else {
        f3.denominator = f1.denominator * f2.denominator;
        f3.numerator = f1.numerator * f2.denominator + f2.numerator * f1.denominator;
    }
    f3.reduce_fraction();
    return f3;
}

fraction operator + (const fraction & f1, int int1) {
    fraction f3;

    if (f1.denominator == 1) {
        f3.numerator = f1.numerator + int1;
    } else {
        f3.numerator = f1.numerator + int1 * f1.denominator;
    }
    return (f3);
}

void operator += (fraction & f1,
    const fraction & f2) {
    f1 = f1 + f2;
}

fraction operator - (const fraction & f1,
    const fraction & f2) {
    fraction f3;
    if (f1.denominator == f2.denominator) {
        f3.numerator = f1.numerator - f2.numerator;
        f3.denominator = f1.denominator;
    } else {
        f3.denominator = f1.denominator * f2.denominator;
        f3.numerator = f1.numerator * f2.denominator - f2.numerator * f1.denominator;
        f3.reduce_fraction();
    }
    return f3;
}

fraction operator - (const fraction & f1, int int1) {
    fraction f2(int1);
    fraction f3;
    f3 = f1 - f2; //use the fraction-fraction overloaded operator
    return (f3);
}

void operator -= (fraction & f1,
    const fraction & f2) {
    f1 = f1 - f2;
}

fraction operator * (const fraction & f1,
    const fraction & f2) {
    fraction f3;
    f3.denominator = f1.denominator * f2.denominator;
    f3.numerator = f1.numerator * f2.numerator;
    f3.reduce_fraction();
    return f3;
}

fraction operator * (const fraction & f1, int int1) {
    fraction f2(int1);
    fraction f3;
    f3 = f1 * f2; //use the fraction*fraction overloaded operator
    return (f3);
}

void operator *= (fraction & f1,
    const fraction & f2) {
    f1 = f1 * f2;
}

fraction operator / (const fraction & f1,
    const fraction & f2) {
    if (f2.numerator == 0) {
        exit(3);
    }
    fraction f3;
    f3.denominator = f1.denominator * f2.numerator;
    f3.numerator = f1.numerator * f2.denominator;
    f3.reduce_fraction();
    return f3;
}

fraction operator / (const fraction & f1, int int1) {
    fraction f2(int1);
    fraction f3;
    f3 = f1 / f2; //use the fraction+fraction overloaded operator
    return (f3);
}
void operator /= (fraction & f1,
    const fraction & f2) {
    f1 = f1 / f2;
}

std::ostream & operator << (std::ostream & os,
    const fraction & f1) {
    if (f1.denominator == 1 || f1.numerator == 0) {
        os << f1.numerator;
        return (os);
    } else if (f1.numerator == f1.denominator) {
        os << 1;
        return (os);
    } else if (f1.denominator == 0) {
        os << "CANNOT DIVIDE BY 0";
        return (os);
    } else {
        std::string fraction_string = to_string(f1.numerator) + '/' + to_string(f1.denominator);
        os << fraction_string;
    }
    return os;
}
bool operator!(const fraction & f) {
    //denominator doesn't matter, only the numeratorerator. 
    return (!(f.numerator));
}
bool operator != (const fraction & f1,
    const fraction & f2) {
    //well the fraction would have already been simplified, so. 
    //so we check to see if both the numerators and denominators are equal to eachother.
    return (!((f1.numerator == f2.numerator) && (f1.denominator == f2.denominator)));
}
bool operator != (const fraction & f1, int int1) {
    //wait so we only need to check if the numerator==int 1 and if f1's denominator is 1, since f1 should have already been simplified.
    return (!((f1.numerator == int1) && (f1.denominator == 1)));
}
bool operator == (const fraction & f1,
    const fraction & f2) {
    return ((f1.numerator == f2.numerator) && (f1.denominator == f2.denominator));
}
bool operator == (const fraction & f1, int int1) {
    return (((f1.numerator == int1) && (f1.denominator == 1)));
}

//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */

class vector2D {

    public:
    vector2D();
    vector2D(int, int);
    vector2D(int, int,  std::vector < fraction > elements);
    vector2D(const std::vector <std::vector <fraction> > &vector_of_vectors );

    bool is_empty();

    int get_number_of_rows() const;
    int get_number_of_columns() const;

    void set_number_of_rows(int);
    void set_number_of_columns(int);

    void keep_rows(const std::vector < int >& rows_to_keep);
    void keep_rows(int);
    void keep_columns(const std::vector < int >& columns_to_keep);
    void flip();
    void row_swap(int, int);
    void column_swap(int, int);
    void multiply_row(int, const fraction&);
    void divide_row(int, const fraction&);
    void multiply_column(int, const fraction&);
    void divide_column(int, const fraction&);

    fraction get_element(int, int) const;

    std::vector < std::vector <fraction> > get_table() const;
    void set_table(const std::vector < std::vector < fraction > > & new_table);
    void set_element(int, int, const fraction&);

    friend std::ostream & operator << (std::ostream & os, const vector2D & vector_2D);
    void print_table();

    protected: int number_of_rows;
    int number_of_columns;

    private: 
    std::vector < std::vector < fraction > > table;

};

    vector2D::vector2D()
        : number_of_rows(0),number_of_columns(0),table(std::vector<std::vector<fraction> >()) {

        }
        
        vector2D::vector2D(int rows, int columns, std::vector < fraction > elements)
    : number_of_rows(rows), number_of_columns(columns) {
            table.resize(number_of_rows, std::vector<fraction>(number_of_columns)); // Initialize the 2D vector
            for (int i = 0; i < number_of_rows; i++) {
                for (int j = 0; j < number_of_columns; j++) {
                    table[i][j] = elements.back();
                    elements.pop_back();
                }
            }
}

  void vector2D::set_number_of_rows(int rows){
    number_of_rows=rows;
  }

    void vector2D:: set_number_of_columns(int cols){
       number_of_columns=cols;
    }

void vector2D:: set_element(int row, int col, const fraction& new_value){
 if (row >= 0 && row < number_of_rows && col >= 0 && col < number_of_columns) {
        table[row][col] = new_value;
    } else {
        exit(6);
    }
}


vector2D :: vector2D(const std:: vector <  std::vector <fraction> >& vector_of_vectors)
    : number_of_rows(vector_of_vectors.size()), number_of_columns(vector_of_vectors[0].size()),table (vector_of_vectors){

}

vector2D:: vector2D(int rows, int cols)
    : number_of_rows(rows), number_of_columns(cols) {
    table.resize(rows, std::vector<fraction>(cols)); // Initialize the 2D vector with default fractions
    }

bool vector2D::is_empty() {
    if (table.empty()) {
        return true;
    }
    return false;
}

int vector2D::get_number_of_rows() const {
    return number_of_rows;
}

int vector2D::get_number_of_columns()const {
    return number_of_columns;
}


void vector2D::keep_columns(const std::vector < int > &columns_to_keep) {
    flip();
    int columns_kept = columns_to_keep.size();
    int j = 0;
    bool keep_column = false;

    while (j < number_of_columns) {
        for (int k = 0; k < columns_kept; k++) {
            if (j == columns_to_keep[k]) {
                keep_column = true;
                break;
            }
        }
        if (!keep_column) {
            keep_column = false;
            table.erase(table.begin() + j);
        } else {
            j++;
        }
    }
    flip();

    number_of_columns= columns_kept;
}

void vector2D::keep_rows(const std::vector < int > &rows_to_keep) {

    int rows_kept = rows_to_keep.size();
    int i = 0;
    bool keep_row = false;
    while (i < number_of_rows) {
        for (int k = 0; k < rows_kept; k++) {
            if (i == rows_to_keep[k]) {
                keep_row = true;
                break;
            }
        }

        if (!keep_row) {
            keep_row = false;
            table.erase(table.begin() + i);
        } else {
            i++;
        }
    }
    number_of_rows=rows_kept;
}

void vector2D::flip() {
    vector2D flipped_vector2D(number_of_columns, number_of_rows);
    //flip 2d vector to column-major order
    for (int i = 0; i < number_of_rows; i++) {
        for (int j = 0; j < number_of_columns; j++) {
            flipped_vector2D.table[j][i] = table[i][j];
        }
    }
    table = flipped_vector2D.table;
    int temp = number_of_columns;
    number_of_columns = number_of_rows;
    number_of_rows = temp;
}

void vector2D::row_swap(int a, int b) {
    if (a >= number_of_rows || b >= number_of_rows) {
        exit(1);
    } else if (a < 0 || b < 0) {
        exit(2);
    } else if (a == b) {
        return;
    }

    std::vector < fraction > temp;

    temp = table[a];
    table[a] = table[b];
    table[b] = temp;
    return;
}

void vector2D::column_swap(int a, int b) {
    if (a >= number_of_columns || b >= number_of_columns) {
        exit(1);
    } else if (a < 0 || b < 0) {
        exit(2);
    } else if (a == b) {
        return;
    }

    for (int i = 0; i < number_of_rows; i++) {
        fraction temp = table[i][a];
        table[i][a] = table[i][b];
        table[i][b] = temp;
    }

}

void vector2D::multiply_row(int a, const fraction &c) {
    for (int j = 0; j < number_of_columns; j++) {
        if (table[a][j] == 0) {
            continue;
        }
        table[a][j] *= c;
    }
    return;
}

void vector2D::multiply_column(int a, const fraction &c) {
    for (int i = 0; i < number_of_rows; i++) {
        if (table[i][a].get_numerator()) {
            continue;
        }
        table[i][a] *= c;
    }
    return;
}

void vector2D::divide_row(int a, const fraction &c ) {
    if (!c.get_numerator()) {
        exit(3);
    }

    for (int j = 0; j < number_of_columns; j++) {
        if (table[a][j].get_numerator()) {
            continue;
        }
        table[a][j] /= c;
    }
    return;
}

void vector2D::divide_column(int a, const fraction &c) {
    if (!c.get_numerator()) {
        exit(3);
    }

    for (int i = 0; i < number_of_rows; i++) {
        if (table[i][a].get_numerator()) {
            continue;
        }
        table[i][a] /= c;
    }
    return;
}

fraction vector2D::get_element(int n, int m) const {
    return (table[n][m]);
}

std::vector < std::vector < fraction > > vector2D::get_table() const{
    return (table);
}

void vector2D::set_table(const std::vector < std::vector < fraction > > &new_table) {
    table = new_table;
    number_of_rows= new_table.size();
    number_of_columns= new_table[0].size();
}

void vector2D::keep_rows(int how_many_rows) {

    if (how_many_rows<0)
        exit(4);

    for (int i = number_of_rows; i < how_many_rows; i--) {
        table.pop_back();
    }
    number_of_rows= how_many_rows;
}

std::ostream & operator << (std::ostream & os,const vector2D & vector_2D) {

    int col1 = 10;
    os << "┌";
    os << "┐" << std::endl;
    for (int i = 0; i < vector_2D.number_of_rows; i++) {
        os << "|";
        for (int j = 0; j < vector_2D.number_of_columns; j++) {
            os << std::setw(col1) << std::left << vector_2D.table[i][j];
        }
        os << "|";
        os << std::endl;
    }
    os << "└";
    for (int k = 0; k < vector_2D.number_of_columns * col1; k++) {
        os << " ";
    }
    os << "┘" << std::endl;
    os << std::endl;

    return (os);
    //┌ ┐ └ ┘
}

void vector2D::print_table() {
    int col1 = 10;
    std::cout << "┌";
    for (int k = 0; k < number_of_columns * col1; k++) {
        std::cout << " ";
    }
    std::cout << "┐" << std::endl;
    for (int i = 0; i < number_of_rows; i++) {
        std::cout << "|";
        for (int j = 0; j < number_of_columns; j++) {
            std::cout << std::setw(col1) << std::left << table[i][j];
        }
        std::cout << "|";
        std::cout << std::endl;
    }
    std::cout << "└";
    for (int k = 0; k < number_of_columns * col1; k++) {
        std::cout << " ";
    }
    std::cout << "┘" << std::endl;
    std::cout << std::endl;

}

//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */
//********************************************************************************************************* */

class matrix: public vector2D {
    public:

    //non-augmented
    matrix(int, int, std::vector < fraction > );
    matrix(const std::vector < std::vector < fraction > > & );
    matrix(const vector2D &);

    //augmented
    matrix(int, int, std::vector < fraction > , std::vector < fraction > );
    matrix( const std::vector < std::vector < fraction > > & , const std::vector < std::vector < fraction > > &);
    matrix(const vector2D &, const vector2D&);

    void print_columnspace();
    void print_rowspace();
    void print_reduced_form();

    std::string string_columnspace() ;
    std::string string_rowspace()  ;
    std::string string_reduced_form() ;

    private: 
    vector2D columnspace;
    vector2D rowspace;
    vector2D reduced_form;
    vector2D vectors;

    void row_reduce();
    std::string solution;
    bool augmented;

};

//NON AUGMENTED 1
matrix::matrix(int n, int m,  std::vector < fraction > elements)
    :vector2D(n,m,elements){
        augmented=false;

    }

//NON AUGMENTED 2

matrix::matrix(const std::vector < std::vector < fraction > > & vector_of_vectors)
    :vector2D(vector_of_vectors){
        augmented=false;

}

//NON AUGMENTED 3

matrix::matrix(const vector2D & vector_2D)
    :vector2D(vector_2D.get_table()) {
        augmented=false;
    
}

//AUGMENTED  1  

matrix::matrix(int rows, int cols, std::vector<fraction> elements, std::vector<fraction> vector_elements)
    : vector2D(rows, cols, elements) {
    
    if (vector_elements.size() % rows != 0) {
        exit(5);    }
    
    vectors = vector2D(rows, vector_elements.size() / rows, vector_elements); // Initialize vectors after the check
    augmented=true;
}

    
//AUGMENTED 2 
matrix::matrix(const std::vector < std::vector <fraction > >& vector_of_vectors , const std::vector < std::vector < fraction > >& user_vectors)
    : vector2D(vector_of_vectors), vectors(user_vectors)
    {
        augmented=true;
    }

// AUGMENTED 3
matrix::matrix(const vector2D &the_matrix, const vector2D &the_vectors)
    :vector2D(the_matrix.get_table()), vectors(the_vectors){

        augmented=true;
    }
    


void matrix::row_reduce() {

        if (!solution.empty()) {
            return;
        }

        reduced_form.set_table(get_table()); //set the reduced matrix to the original matrix from the parent.
        std::vector < int > pivot_positions; //stack storing the position of each pivot. The index corresponds to the row, and the value corresponds to the column. Since a row cannot be skipped in row reducing, the pivot rows will always be sequential starting from 0.
        int current_row = 0;
        int current_col = 0;

        while (current_row < number_of_rows && current_col < number_of_columns) {

            //1: ROW SWAP IF CURRENT POSITION IS ZERO
            if (!reduced_form.get_element(current_row, current_col).get_numerator()) {
                bool no_swap = true;

                //checks the next elements in the column until a non zero element is found. Once found, the rows are swapped.
                for (int i = current_row + 1; i < number_of_rows; i++) {
                    if (reduced_form.get_element(i, current_col).get_numerator()) {
                        reduced_form.row_swap(current_row, i);
                        if (augmented)
                            vectors.row_swap(current_row, i);
                        solution += "R" + to_string(current_row) + "<---> R" + to_string(i) + "\n";
                        no_swap = false;
                        break;
                    }
                }
                //Move to the next column
                if (no_swap) {
                    current_col++;
                    continue;
                }
            }

            pivot_positions.push_back(current_col);

            //2: ROW DIVISION IF CURRENT POSITION IS NOT 1
            fraction division_constant = reduced_form.get_element(current_row, current_col);

            if (division_constant != 1) {
                reduced_form.divide_row(current_row, division_constant);
                if (augmented) {
                    vectors.divide_row(current_row, division_constant);
                }
                solution += "R" + to_string(current_row) + "/" + to_string(division_constant) + "\n";
            }

            //3: ELIMINATE ELEMENTS IN THE COLLUMN BELOW THE PIVOT 
            for (int i = current_row + 1; i < number_of_rows; i++) {
                if (!reduced_form.get_element(i,current_col)){
                    continue; //no need to do a row operation
                }
                //always divided by 1 tho fix it
                fraction elimination_constant = (reduced_form.get_element(i, current_col) / reduced_form.get_element(current_row, current_col));
              
                if (elimination_constant.get_numerator() > 0) {
                    solution += "R" + to_string(i) + "-" + to_string(elimination_constant) + "R" + to_string(current_row) + "\n";
                }

                //ensures there are no double negatives in the output
                if (elimination_constant.get_numerator() < 0) {
                    solution += "R" + to_string(i) + "+" + to_string(elimination_constant * -1) + "R" + to_string(current_row) + "\n";
                }

                //adds/subtracts a constant multiple of the current row to each row below to eliminate 
                for (int j = current_col; j < number_of_columns; j++) {
                    fraction new_value = reduced_form.get_element(i, j) - elimination_constant * reduced_form.get_element(current_row, j);
                    reduced_form.set_element(i, j, new_value);
                }

                if (augmented) {
                    for (int j = 0; j < vectors.get_number_of_columns(); j++) {
                        fraction new_value = vectors.get_element(i, j) - elimination_constant * vectors.get_element(current_row, j);
                        vectors.set_element(i, j, new_value);
                    }

                }

            }
            solution += "The matrix is now in Echelon Form. \n";

            // The matrix now needs to be reduced on the upper triangle. In each pivot column, The elements above a pivot
            int pivots_left = pivot_positions.size()-1;

            while (pivots_left >= 0) {
                current_col = pivot_positions.back();
                current_row = pivots_left;

                for (int i = current_row - 1; i >= 0; i--) {
                    if (!reduced_form.get_element(i,current_col)){
                    continue; //no need to do a row operation
                }
                    fraction elimination_constant = reduced_form.get_element(i, current_col) / reduced_form.get_element(current_row, current_col);

                    if (elimination_constant.get_numerator() > 0) {
                        solution += "R" + to_string(i) + "-" + to_string(elimination_constant) + "R" + to_string(current_row) + "\n";
                    }

                    if (elimination_constant.get_numerator() < 0) {
                        solution += "R" + to_string(i) + "+" + to_string(elimination_constant * -1) + "R" + to_string(current_row) + "\n";
                    }

                    for (int j = current_col; j >= 0; j++) {
                        fraction new_value = reduced_form.get_element(i, j) - elimination_constant * reduced_form.get_element(current_row, j);
                        reduced_form.set_element(i, j, new_value);
                    }
                    if (augmented) {
                        for (int j = 0; j < vectors.get_number_of_columns(); j++) {
                            fraction new_value = vectors.get_element(i, j) - elimination_constant * vectors.get_element(current_row, j);
                            vectors.set_element(i, j, new_value);

                        }
                    }

                }
                pivots_left--;
            }
        }

            //Extracting the results
            solution += "The Matrix is in Reduced Row Echelon Form. \n";
            solution += "-> Number of indendant columns: " + to_string(pivot_positions.size()) + "\n";
            solution += "-> Number of dependant columns: " + to_string(number_of_columns - pivot_positions.size()) + "\n";

            //The columnspace is the set of independant columns from the original matrix
            columnspace.set_table(get_table());
            columnspace.keep_columns(pivot_positions);

            //The row space is the set of independant rows from the reduced matrix
            rowspace.set_table(reduced_form.get_table());
            rowspace.keep_rows(pivot_positions.size());

            solution += "-> ROWSPACE: \n";
            solution += string_rowspace();

            solution += solution += "-> COLLUMNSPACE: \n";
            solution += string_columnspace();
        }

        void matrix::print_rowspace() {

            if (solution.empty()) {
                row_reduce();
            }

            int rowspace_rows = rowspace.get_number_of_rows();
            int rowspace_columns = rowspace.get_number_of_columns();

            for (int i = 0; i < rowspace_rows; i++) {
                std::cout << "{ ";
                for (int j = 0; j < rowspace_columns; j++) {
                    std::cout << reduced_form.get_element(i, j);
                    if (i != rowspace_rows - 1)
                        std::cout << ", ";
                }
                std::cout << "}" << std::endl;
            }

        }

        void matrix::print_reduced_form() {
            if (solution.empty()) {
                row_reduce();
            }

            if (augmented) {
                reduced_form.print_table();
            } else {
                int col1 = 10;
                std::cout << "┌";
                for (int k = 0; k < (number_of_columns + vectors.get_number_of_columns()) * col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┐" << std::endl;
                for (int i = 0; i < number_of_rows; i++) {
                    std::cout << "|";
                    for (int j = 0; j < number_of_columns; j++) {
                        std::cout << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    }
                    std::cout << "|";
                    for (int j = 0; j < vectors.get_number_of_columns(); j++) {
                        std::cout << std::setw(col1) << std::left << vectors.get_element(i, j);
                        std::cout << "|";
                    }
                    std::cout << std::endl;
                }
                std::cout << "└";
                for (int k = 0; k < (number_of_columns + vectors.get_number_of_columns()) * col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┘" << std::endl;
                std::cout << std::endl;
            }
        }

        void matrix::print_columnspace() {
            if (solution.empty()) {
                row_reduce();
            }

            int col1 = 10;
            for (int i = 0; i < columnspace.get_number_of_columns(); i++) {
                std::cout << "┌";
                for (int k = 0; k < col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┐";
                std::cout << "  ";
            }

            for (int i = 0; i < columnspace.get_number_of_columns(); i++) {
                std::cout << "|";
                for (int j = 0; j < number_of_columns; j++) {
                    std::cout << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    std::cout << "|";
                    std::cout << "  ";
                }
                std::cout << std::endl;
            }
            for (int i = 0; i < columnspace.get_number_of_columns(); i++) {
                std::cout << "└";
                for (int k = 0; k < col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┘";
            }

            std::cout << std::endl;

        }

        std::string matrix::string_columnspace() {
            if (solution.empty()) {
                row_reduce();
            }

            std::ostringstream os;
            int col1 = 10;
            for (int i = 0; i < columnspace.get_number_of_columns(); i++) {
                os << "┌";
                for (int k = 0; k < col1; k++) {
                    os << " ";
                }
                os << "┐";
                os << "  ";
            }

            for (int i = 0; i < columnspace.get_number_of_columns(); i++) {
                os << "|";
                for (int j = 0; j < number_of_columns; j++) {
                    os << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    os << "|";
                    os << "  ";
                }
                os << std::endl;
            }
            for (int i = 0; i < columnspace.get_number_of_columns(); i++) {
                os << "└";
                for (int k = 0; k < col1; k++) {
                    os << " ";
                }
                os << "┘";
            }

            os << std::endl;
            return os.str();

        }

        std::string matrix::string_rowspace()  {
            if (solution.empty()) {
                row_reduce();
            }

            std::ostringstream oss;
            int rowspace_rows = rowspace.get_number_of_rows();
            int rowspace_columns = rowspace.get_number_of_columns();

            for (int i = 0; i < rowspace_rows; i++) {
                oss << "{ ";
                for (int j = 0; j < rowspace_columns; j++) {
                    oss << reduced_form.get_element(i, j);
                    if (j != rowspace_columns - 1)
                        oss << ", ";
                }
                oss << "}" << std::endl;
            }

            return oss.str();
        }

        std::string matrix::string_reduced_form() {
            if (solution.empty()) {
                row_reduce();
            }

            std::ostringstream oss;

            if (vectors.is_empty()) {
                int col1 = 10;
                oss << "┌";
                for (int k = 0; k < number_of_columns * col1; k++) {
                    oss << " ";
                }
                oss << "┐" << std::endl;
                for (int i = 0; i < number_of_rows; i++) {
                    oss << "|";
                    for (int j = 0; j < number_of_columns; j++) {
                        oss << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    }
                    oss << "|";
                    oss << std::endl;
                }
                oss << "└";
                for (int k = 0; k < number_of_columns * col1; k++) {
                    oss << " ";
                }
                oss << "┘" << std::endl;
                oss << std::endl;
                return oss.str();
            } else {
                int col1 = 10;
                oss << "┌";
                for (int k = 0; k < (number_of_columns + vectors.get_number_of_columns()) * col1; k++) {
                    oss << " ";
                }
                oss << "┐" << std::endl;

                for (int i = 0; i < number_of_rows; i++) {
                    oss << "|";
                    for (int j = 0; j < number_of_columns; j++) {
                        oss << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    }
                    oss << "|";
                    for (int j = 0; j < vectors.get_number_of_columns(); j++) {
                        oss << std::setw(col1) << std::left << vectors.get_element(i, j);
                        oss << "|";
                    }
                    oss << std::endl;
                }

                oss << "└";
                for (int k = 0; k < (number_of_columns + vectors.get_number_of_columns()) * col1; k++) {
                    oss << " ";
                }
                oss << "┘" << std::endl;
                oss << std::endl;
            }

            return oss.str();
        }

    int main(){
        std::cout<<"hello";
        return 0;
    }
