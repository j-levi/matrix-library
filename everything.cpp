#include <vector>

#include <iomanip>

#include <sstream>

#include <iostream>

#include <string>

#include <algorithm>

#include <deque>

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

class fraction {
    public: int get_numerator();
    fraction(); // sets numerator to 0 and denominator to 1
    fraction(int); // sets numerator to the paramamter and denominator to 1
    fraction(int, int); //sets the numerator and denominator to the paramaters
    fraction(const std::string & input); // sets the fraction to a string
    void reduce_fraction(); // reduces fraction to its simplest form.The function divides the numerator and denomitaor by the GCD.
    int gcd_Euclidean(int, int); //helper function for the reduce_fraction. determinbes the Greatest Common Divisor of two integers.
    friend int length(const fraction & ); //finds the length of the fraction using string stream; useful for printing the fraction in a matrix to calculate padding.
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

    friend std::vector <fraction> input_vector(int n, int m);

    private: int numerator;
    int denominator;
};

std::vector <fraction>  input_vector(int n, int m){

        int col1=10;

        std::cout<< "┌";
        print_spaces(m* col1);
        std::cout<<"┐";
        std::cout<<std::endl;
        std::vector <fraction> vector_input; 

    
        for (int i=0; i<n;i++)
        {
            std::cout <<"|";
            for (int j=0; j<m; j++)
            {
                fraction input;
                std::cin>> input;
                vector_input.push_back(input); 
                int element_length= length(input);
                print_spaces(col1-element_length);
            }
            std::cout <<"|";
            std::cout<<std::endl;
        }
        return(vector_input);
}


int fraction::get_numerator() {
    return numerator;
}

fraction::fraction() {
    numerator = 0;
    denominator = 1;
}

fraction::fraction(int n): numerator(n), denominator(1) {}

fraction::fraction(int n, int d): numerator(n), denominator(d) {
    if (denominator == 0) {
        throw std::invalid_argument("Denominator cannot be zero");
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
                throw std::invalid_argument("Denominator cannot be zero");
            }
        } else {
            denominator = 1;
        }
    } else {
        throw std::invalid_argument("Invalid input format");
    }
    reduce_fraction();
}

int fraction::gcd_Euclidean(int a, int b) {
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
    std::string p1 = std::to_string(f.numerator);
    std::string p2 = std::to_string(f.denominator);
    return (p1.length() + 1 + p2.length());
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
        throw std::invalid_argument("Invalid input format");
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
        throw std::invalid_argument("Division by zero");
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

class vector2D {

    public:
    vector2D();

    vector2D(int, int);
    vector2D(int, int, std::vector < fraction > elements);
    vector2D(int, int, std::vector <std::vector <fraction> >vector_of_vectors );

    bool is_empty();

    int get_number_of_rows();
    int get_number_of_columns();

    void keep_rows(std::vector < int > rows_to_keep);
    void keep_rows(int);
    void keep_columns(std::vector < int > collumns_to_keep);
    void flip();
    void row_swap(int, int);
    void collumn_swap(int, int);
    void multiply_row(int, fraction);
    void divide_row(int, fraction);
    void multiply_column(int, fraction);
    void divide_column(int, fraction);

    void multiply_column_below(int, int, fraction);

    fraction get_element(int, int);

    std::vector < std::vector <fraction> > get_table();
    void set_table(std::vector < std::vector < fraction > > new_table);
    void set_element(int, int, fraction);

    friend std::ostream & operator << (std::ostream & os,
    const vector2D & vector_2D);
    void print_table();

    protected: int number_of_rows;
    int number_of_collumns;

    private: std::vector < std::vector < fraction > > table;

};

vector2D::vector2D() {

}

vector2D::vector2D(int rows, int collumns, std::vector < fraction > elements) {
    number_of_rows = rows;
    number_of_collumns = collumns;
    for (int i = 0; i < number_of_rows; i++) {
        for (int j = 0; j < number_of_collumns; j++) {
            table[i][j] = elements.back();
            elements.pop_back();
        }
    }
}

vector2D :: vector2D(int n, int m, std:: vector < std::vector <fraction> > vector_of_vectors){
    number_of_rows=n;
    number_of_collumns=m;
    table = vector_of_vectors;
}

bool vector2D::is_empty() {
    if (table.empty()) {
        return true;
    }
    return false;
}

int vector2D::get_number_of_rows() {
    return number_of_rows;
}

int vector2D::get_number_of_columns() {
    return number_of_collumns;
}


void vector2D::keep_columns(std::vector < int > collumns_to_keep) {
    flip();
    int collumns_kept = collumns_to_keep.size();
    int j = 0;
    bool keep_collumn = false;

    while (j < number_of_collumns) {
        for (int k = 0; k < collumns_kept; k++) {
            if (j == collumns_to_keep[k]) {
                keep_collumn = true;
                break;
            }
        }
        if (!keep_collumn) {
            keep_collumn = false;
            table.erase(table.begin() + j);
        } else {
            j++;
        }
    }
    flip();
}

void vector2D::keep_rows(std::vector < int > rows_to_keep) {

    int rows_kept = rows_to_keep.size();
    int i = 0;
    bool keep_row = false;
    while (i < number_of_collumns) {
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
}

void vector2D::flip() {
    vector2D flipped_vector2D(number_of_collumns, number_of_rows);
    //flip 2d vector to column-major order
    for (int i = 0; i < number_of_rows; i++) {
        for (int j = 0; j < number_of_collumns; j++) {
            flipped_vector2D.table[j][i] = table[i][j];
        }
    }
    table = flipped_vector2D.table;
    int temp = number_of_collumns;
    number_of_collumns = number_of_rows;
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

void vector2D::collumn_swap(int a, int b) {
    if (a >= number_of_collumns || b >= number_of_collumns) {
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

void vector2D::multiply_row(int a, fraction c) {
    for (int j = 0; j < number_of_collumns; j++) {
        if (table[a][j] == 0) {
            continue;
        }
        table[a][j] *= c;
    }
    return;
}

void vector2D::multiply_column(int a, fraction c) {
    for (int i = 0; i < number_of_rows; i++) {
        if (table[i][a].get_numerator()) {
            continue;
        }
        table[i][a] *= c;
    }
    return;
}

void vector2D::divide_row(int a, fraction c) {
    if (c.get_numerator()) {
        exit(3);
    }

    for (int j = 0; j < number_of_collumns; j++) {
        if (table[a][j].get_numerator()) {
            continue;
        }
        table[a][j] /= c;
    }
    return;
}

void vector2D::divide_column(int a, fraction c) {
    if (c.get_numerator()) {
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

fraction vector2D::get_element(int n, int m) {
    return (table[n][m]);
}

std::vector < std::vector < fraction > > vector2D::get_table() {
    return (table);
}

void vector2D::set_table(std::vector < std::vector < fraction > > new_table) {
    table = new_table;

}

void vector2D::multiply_column_below(int col, int row, fraction constant) {
    for (int i = row + 1; i < number_of_rows; i++) {
        table[i][col] *= constant;
    }
}

void vector2D::keep_rows(int how_many_rows) {

    if (how_many_rows<0)
        exit(4);

    for (int i = number_of_rows; i < how_many_rows; i--) {
        table.pop_back();
    }
}

std::ostream & operator << (std::ostream & os,const vector2D & vector_2D) {

    int col1 = 10;
    os << "┌";
    os << "┐" << std::endl;
    for (int i = 0; i < vector_2D.number_of_rows; i++) {
        os << "|";
        for (int j = 0; j < vector_2D.number_of_collumns; j++) {
            os << std::setw(col1) << std::left << vector_2D.table[i][j];
        }
        os << "|";
        os << std::endl;
    }
    os << "└";
    for (int k = 0; k < vector_2D.number_of_collumns * col1; k++) {
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
    for (int k = 0; k < number_of_collumns * col1; k++) {
        std::cout << " ";
    }
    std::cout << "┐" << std::endl;
    for (int i = 0; i < number_of_rows; i++) {
        std::cout << "|";
        for (int j = 0; j < number_of_collumns; j++) {
            std::cout << std::setw(col1) << std::left << table[i][j];
        }
        std::cout << "|";
        std::cout << std::endl;
    }
    std::cout << "└";
    for (int k = 0; k < number_of_collumns * col1; k++) {
        std::cout << " ";
    }
    std::cout << "┘" << std::endl;
    std::cout << std::endl;

}


class matrix: public vector2D {
    public:

    //non-augmented
    matrix(int, int, std::vector < fraction > );
    matrix(int, int, std::vector < std::vector < fraction > > );
    matrix(vector2D);

    //augmented
    matrix(int, int, std::vector < fraction > , std::vector < fraction > );
    matrix(int, int, std::vector < std::vector < fraction > > , std::vector < std::vector < fraction > > );
    matrix(vector2D, vector2D);

    void print_collumnspace();
    void print_rowspace();
    void print_reduced_form();

    std::string string_collumnspace();
    std::string string_rowspace();
    std::string string_reduced_form();

    private: vector2D collumnspace;
    vector2D rowspace;
    vector2D reduced_form;
    vector2D vectors;

    void row_reduce();
    std::string solution;

};

//NON AUGMENTED
matrix::matrix(int n, int m, std::vector < fraction > elements) {
    vector2D(n,m,elements);
}

matrix::matrix(int rows, int cols, std::vector < std::vector < fraction > > vector_of_vectors){
    vector2D(rows, cols, vector_of_vectors);

}

matrix::matrix(vector2D vector_2D) {
    number_of_rows = vector_2D.get_number_of_rows();
    number_of_collumns = vector_2D.get_number_of_columns();
    set_table(vector_2D.get_table());
}

//AUGMENTED
matrix::matrix(vector2D the_matrix, vector2D the_vectors) {

    number_of_rows = the_matrix.get_number_of_rows();
    number_of_collumns = the_matrix.get_number_of_columns();

    set_table(the_matrix.get_table());
    the_vectors.flip();
    vectors.set_table(the_vectors.get_table());

}

void matrix::row_reduce() {

        bool augmented = false;

        if (vectors.is_empty()) {
            augmented = true;
        }

        if (!solution.empty()) {
            return;
        }

        reduced_form.set_table(get_table()); //set the reduced matrix to the original matrix from the parent.
        std::vector < int > pivot_positions; //stack storing the position of each pivot. The index corresponds to the row, and the value corresponds to the column. Since a row cannot be skipped in row reducing, the pivot rows will always be sequential starting from 0.
        int current_row = 0;
        int current_col = 0;

        while (current_row < number_of_rows && current_col < number_of_collumns) {

            //1: ROW SWAP IF CURRENT POSITION IS ZERO
            if (!reduced_form.get_element(current_row, current_col).get_numerator()) {
                bool no_swap = false;

                //checks the next elements in the collumn until a non zero element is found. Once found, the rows are swapped.
                for (int i = current_row + 1; i < number_of_rows; i++) {
                    if (reduced_form.get_element(i, current_col).get_numerator()) {
                        reduced_form.row_swap(current_row, i);
                        if (augmented)
                            vectors.row_swap(current_row, i);
                        solution += "R" + to_string(current_row) + "<---> R" + to_string(i) + "/n";
                        no_swap = true;
                        break;
                    }
                }
                //Move to the next collumn
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
                solution += "R" + to_string(current_row) + "/" + to_string(division_constant) + "/n";
            }

            //3: ELIMINATE ELEMENTS IN THE COLLUMN BELOW THE PIVOT 
            for (int i = current_row + 1; i < number_of_rows; i++) {
                fraction elimination_constant = (reduced_form.get_element(i, current_col) / reduced_form.get_element(current_row, current_col));

                if (elimination_constant.get_numerator() > 0) {
                    solution += "R" + to_string(i) + "-" + to_string(elimination_constant) + "R" + to_string(current_row) + "/n";
                }

                //ensures there are no double negatives in the output
                if (elimination_constant.get_numerator() < 0) {
                    solution += "R" + to_string(i) + "+" + to_string(elimination_constant * -1) + "R" + to_string(current_row) + "/n";
                }

                //adds/subtracts a constant multiple of the current row to each row below to eliminate 
                for (int j = current_col; j < number_of_collumns; j++) {
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
            solution += "The matrix is now in Echelon Form. /n";

            // The matrix now needs to be reduced on the upper triangle. In each pivot collumn, The elements above a pivot
            int pivots_left = pivot_positions.size();

            while (pivots_left >= 0) {
                current_col = pivot_positions.back();
                current_row = pivots_left;

                for (int i = current_row - 1; i >= 0; i--) {
                    fraction elimination_constant = reduced_form.get_element(i, current_col) / reduced_form.get_element(current_row, current_col);

                    if (elimination_constant.get_numerator() > 0) {
                        solution += "R" + to_string(i) + "-" + to_string(elimination_constant) + "R" + to_string(current_row) + "/n";
                    }

                    if (elimination_constant.get_numerator() < 0) {
                        solution += "R" + to_string(i) + "+" + to_string(elimination_constant * -1) + "R" + to_string(current_row) + "/n";
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
            solution += "The Matrix is in Reduced Row Echelon Form. /n";
            solution += "-> Number of indendant collumns: " + to_string(pivot_positions.size()) + "/n";
            solution += "-> Number of dependant collumns: " + to_string(number_of_collumns - pivot_positions.size()) + "/n";

            //The collumnspace is the set of independant collumns from the original matrix
            collumnspace.set_table(get_table());
            collumnspace.keep_columns(pivot_positions);

            //The row space is the set of independant rows from the reduced matrix
            rowspace.set_table(reduced_form.get_table());
            rowspace.keep_rows(pivot_positions.size());

        }

        void matrix::print_rowspace() {

            if (solution.empty()) {
                row_reduce();
            }

            int rowspace_rows = rowspace.get_number_of_rows();
            int rowspace_collumns = rowspace.get_number_of_columns();

            for (int i = 0; i < rowspace_rows; i++) {
                std::cout << "{ ";
                for (int j = 0; j < rowspace_collumns; j++) {
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

            if (vectors.is_empty()) {
                reduced_form.print_table();
            } else {
                int col1 = 10;
                std::cout << "┌";
                for (int k = 0; k < (number_of_collumns + vectors.get_number_of_columns()) * col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┐" << std::endl;
                for (int i = 0; i < number_of_rows; i++) {
                    std::cout << "|";
                    for (int j = 0; j < number_of_collumns; j++) {
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
                for (int k = 0; k < (number_of_collumns + vectors.get_number_of_columns()) * col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┘" << std::endl;
                std::cout << std::endl;
            }
        }

        void matrix::print_collumnspace() {
            if (solution.empty()) {
                row_reduce();
            }

            int col1 = 10;
            for (int i = 0; i < collumnspace.get_number_of_columns(); i++) {
                std::cout << "┌";
                for (int k = 0; k < col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┐";
                std::cout << "  ";
            }

            for (int i = 0; i < collumnspace.get_number_of_columns(); i++) {
                std::cout << "|";
                for (int j = 0; j < number_of_collumns; j++) {
                    std::cout << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    std::cout << "|";
                    std::cout << "  ";
                }
                std::cout << std::endl;
            }
            for (int i = 0; i < collumnspace.get_number_of_columns(); i++) {
                std::cout << "└";
                for (int k = 0; k < col1; k++) {
                    std::cout << " ";
                }
                std::cout << "┘";
            }

            std::cout << std::endl;

        }

        std::string matrix::string_collumnspace() {
            if (solution.empty()) {
                row_reduce();
            }

            std::ostringstream os;
            int col1 = 10;
            for (int i = 0; i < collumnspace.get_number_of_columns(); i++) {
                os << "┌";
                for (int k = 0; k < col1; k++) {
                    os << " ";
                }
                os << "┐";
                os << "  ";
            }

            for (int i = 0; i < collumnspace.get_number_of_columns(); i++) {
                os << "|";
                for (int j = 0; j < number_of_collumns; j++) {
                    os << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    os << "|";
                    os << "  ";
                }
                os << std::endl;
            }
            for (int i = 0; i < collumnspace.get_number_of_columns(); i++) {
                os << "└";
                for (int k = 0; k < col1; k++) {
                    os << " ";
                }
                os << "┘";
            }

            os << std::endl;
            return os.str();

        }

        std::string matrix::string_rowspace() {
            if (solution.empty()) {
                row_reduce();
            }

            std::ostringstream oss;
            int rowspace_rows = rowspace.get_number_of_rows();
            int rowspace_collumns = rowspace.get_number_of_columns();

            for (int i = 0; i < rowspace_rows; i++) {
                oss << "{ ";
                for (int j = 0; j < rowspace_collumns; j++) {
                    oss << reduced_form.get_element(i, j);
                    if (j != rowspace_collumns - 1)
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
                for (int k = 0; k < number_of_collumns * col1; k++) {
                    oss << " ";
                }
                oss << "┐" << std::endl;
                for (int i = 0; i < number_of_rows; i++) {
                    oss << "|";
                    for (int j = 0; j < number_of_collumns; j++) {
                        oss << std::setw(col1) << std::left << reduced_form.get_element(i, j);
                    }
                    oss << "|";
                    oss << std::endl;
                }
                oss << "└";
                for (int k = 0; k < number_of_collumns * col1; k++) {
                    oss << " ";
                }
                oss << "┘" << std::endl;
                oss << std::endl;
                return oss.str();
            } else {
                int col1 = 10;
                oss << "┌";
                for (int k = 0; k < (number_of_collumns + vectors.get_number_of_columns()) * col1; k++) {
                    oss << " ";
                }
                oss << "┐" << std::endl;

                for (int i = 0; i < number_of_rows; i++) {
                    oss << "|";
                    for (int j = 0; j < number_of_collumns; j++) {
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
                for (int k = 0; k < (number_of_collumns + vectors.get_number_of_columns()) * col1; k++) {
                    oss << " ";
                }
                oss << "┘" << std::endl;
                oss << std::endl;
            }

            return oss.str();
        }

    int main(){
        std::vector <fraction> input= input_vector(5,5);
        matrix my_matrix(4,4,input);

        return 0;
    }

    
