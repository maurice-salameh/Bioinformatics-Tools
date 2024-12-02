import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Scanner;

public class NeedlemanWunsch {

    public double insertion = -1;
    public double deletion = -1;
    public double mismatch = -1;
    public double match = 1;
    public double gap_open_penalty = -10;
    public double gap_extension_penalty = -0.5;

    public String sequence1;
    public String sequence2;

    public NeedlemanWunsch(String s1, String s2) {
        sequence1 = s1;
        sequence2 = s2;

        Scanner s = new Scanner(System.in);
        try {
            System.out.print("Enter insertion score (default -1): ");
            insertion = getOrDefault(s, -1);

            System.out.print("Enter deletion score (default -1): ");
            deletion = getOrDefault(s, -1);

            System.out.print("Enter mismatch score (default -1): ");
            mismatch = getOrDefault(s, -1);

            System.out.print("Enter match score (default +1): ");
            match = getOrDefault(s, 1);

            System.out.print("Enter gap open penalty (default -10): ");
            gap_open_penalty = getOrDefault(s, -10);

            System.out.print("Enter gap extension penalty (default -0.5): ");
            gap_extension_penalty = getOrDefault(s, -0.5);
        } catch (Exception e) {
            System.out.println("Invalid input. Using default values.");
        }
    }



    //this function does not take into consideration the gap extension and gap opening penalty in the alignment
    public void align() {
        int n = sequence1.length();
        int m = sequence2.length();

        double[][] scoreMatrix = new double[n + 1][m + 1];
        int[][] tracebackMatrix = new int[n + 1][m + 1]; // To trace the alignment path
        scoreMatrix[0][0] = 0;

        // Initialize first column (sequence1 gaps)
        for (int i = 1; i <= n; i++) {
            scoreMatrix[i][0] = scoreMatrix[i - 1][0] + deletion;
            tracebackMatrix[i][0] = 1; // Indicates gap in sequence2
        }

        // Initialize first row (sequence2 gaps)
        for (int j = 1; j <= m; j++) {
            scoreMatrix[0][j] = scoreMatrix[0][j - 1] + insertion;
            tracebackMatrix[0][j] = 2; // Indicates gap in sequence1
        }

        // Fill score matrix
        for (int i = 1; i <= n; i++) {
            char char1 = sequence1.charAt(i - 1);
            for (int j = 1; j <= m; j++) {
                char char2 = sequence2.charAt(j - 1);

                //moving diagonally means there's a nucl at each column and row
                double scoreDiagonal = scoreMatrix[i - 1][j - 1] + ((char1 == char2) ? match : mismatch);

                //moving vertically in the matrix means deletion
                double scoreUpwards = scoreMatrix[i - 1][j] + deletion;

                //moving horizontally in the matrix means insertion
                double scoreLeft = scoreMatrix[i][j - 1] + insertion;

                double maxScore = Math.max(Math.max(scoreDiagonal, scoreUpwards), scoreLeft);
                scoreMatrix[i][j] = maxScore;

                if (maxScore == scoreDiagonal) {
                    tracebackMatrix[i][j] = 0; // Match/Mismatch
                } else if (maxScore == scoreUpwards) {
                    tracebackMatrix[i][j] = 1; // Gap in sequence2
                } else {
                    tracebackMatrix[i][j] = 2; // Gap in sequence1
                }
            }
        }

        // Traceback
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        int i = n, j = m;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && tracebackMatrix[i][j] == 0) {
                alignedSeq1.append(sequence1.charAt(i - 1));
                alignedSeq2.append(sequence2.charAt(j - 1));
                i--;
                j--;
            } else if (i > 0 && tracebackMatrix[i][j] == 1) {
                alignedSeq1.append(sequence1.charAt(i - 1));
                alignedSeq2.append('-');
                i--;
            } else {
                alignedSeq1.append('-');
                alignedSeq2.append(sequence2.charAt(j - 1));
                j--;
            }
        }

        alignedSeq1.reverse();
        alignedSeq2.reverse();

        // Calculate alignment statistics
        int matches = 0, gaps = 0, identities = 0;
        for (int k = 0; k < alignedSeq1.length(); k++) {
            if (alignedSeq1.charAt(k) == alignedSeq2.charAt(k)) {
                identities++;
            }
            if (alignedSeq1.charAt(k) == '-' || alignedSeq2.charAt(k) == '-') {
                gaps++;
            }
        }

        double similarity = (double) identities / alignedSeq1.length() * 100;
        double gapPercent = (double) gaps / alignedSeq1.length() * 100;

        // Print results
        System.out.println("Alignment Result:");
        System.out.println("Aligned Sequence 1: " + alignedSeq1);
        System.out.println("Aligned Sequence 2: " + alignedSeq2);
        System.out.println("Length: " + alignedSeq1.length());
        System.out.printf("Identity: %d/%d (%.2f%%)%n", identities, alignedSeq1.length(), similarity);
        System.out.printf("Gaps: %d/%d (%.2f%%)%n", gaps, alignedSeq1.length(), gapPercent);
        System.out.println("Score: " + scoreMatrix[n][m]);
    }




    private static final double[][] EDNAFULL = {
            { 5, -4, -4, -4 },  // A
            { -4, 5, -4, -4 },  // C
            { -4, -4, 5, -4 },  // G
            { -4, -4, -4, 5 }   // T
    };

    private double getSubstitutionScore(char a, char b) {
        int indexA = "ACGT".indexOf(a);
        int indexB = "ACGT".indexOf(b);
        if (indexA == -1 || indexB == -1) return -1.0; // Non-standard character
        return EDNAFULL[indexA][indexB];
    }
    /*
    public void align_using_affine_gap_penalties() {
        int n = sequence1.length();
        int m = sequence2.length();

        // Initialize matrices
        double[][] M = new double[n + 1][m + 1]; // Match/mismatch matrix
        double[][] X = new double[n + 1][m + 1]; // Gap in sequence2
        double[][] Y = new double[n + 1][m + 1]; // Gap in sequence1
        int[][] tracebackMatrix = new int[n + 1][m + 1]; // To trace the alignment path

        // Initialize the matrices
        M[0][0] = 0;
        for (int i = 1; i <= n; i++) {
            M[i][0] = Double.NEGATIVE_INFINITY;
            X[i][0] = gap_open_penalty + (i - 1) * gap_extension_penalty;
            Y[i][0] = Double.NEGATIVE_INFINITY;
            tracebackMatrix[i][0] = 1; // Gap in sequence2
        }
        for (int j = 1; j <= m; j++) {
            M[0][j] = Double.NEGATIVE_INFINITY;
            X[0][j] = Double.NEGATIVE_INFINITY;
            Y[0][j] = gap_open_penalty + (j - 1) * gap_extension_penalty;
            tracebackMatrix[0][j] = 2; // Gap in sequence1
        }

        // Fill matrices
        for (int i = 1; i <= n; i++) {
            char char1 = sequence1.charAt(i - 1);
            for (int j = 1; j <= m; j++) {
                char char2 = sequence2.charAt(j - 1);

                // Calculate scores for M, X, and Y
                double scoreDiagonal = M[i - 1][j - 1] + getSubstitutionScore(char1, char2);
                M[i][j] = Math.max(scoreDiagonal, Math.max(X[i - 1][j - 1], Y[i - 1][j - 1]));

                double gapInSeq2 = Math.max(M[i - 1][j] + gap_open_penalty, X[i - 1][j] + gap_extension_penalty);
                X[i][j] = gapInSeq2;

                double gapInSeq1 = Math.max(M[i][j - 1] + gap_open_penalty, Y[i][j - 1] + gap_extension_penalty);
                Y[i][j] = gapInSeq1;

                // Update traceback matrix
                if (M[i][j] == scoreDiagonal) {
                    tracebackMatrix[i][j] = 0; // Match/Mismatch
                } else if (M[i][j] == gapInSeq2) {
                    tracebackMatrix[i][j] = 1; // Gap in sequence2
                } else {
                    tracebackMatrix[i][j] = 2; // Gap in sequence1
                }
            }
        }

        // Traceback to build aligned sequences
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        int i = n, j = m;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && tracebackMatrix[i][j] == 0) {
                alignedSeq1.append(sequence1.charAt(i - 1));
                alignedSeq2.append(sequence2.charAt(j - 1));
                i--;
                j--;
            } else if (i > 0 && tracebackMatrix[i][j] == 1) {
                alignedSeq1.append(sequence1.charAt(i - 1));
                alignedSeq2.append('-');
                i--;
            } else {
                alignedSeq1.append('-');
                alignedSeq2.append(sequence2.charAt(j - 1));
                j--;
            }
        }

        alignedSeq1.reverse();
        alignedSeq2.reverse();

        // Calculate alignment statistics
        int identities = 0, gaps = 0;
        for (int k = 0; k < alignedSeq1.length(); k++) {
            if (alignedSeq1.charAt(k) == alignedSeq2.charAt(k)) {
                identities++;
            }
            if (alignedSeq1.charAt(k) == '-' || alignedSeq2.charAt(k) == '-') {
                gaps++;
            }
        }

        double similarity = (double) identities / alignedSeq1.length() * 100;
        double gapPercent = (double) gaps / alignedSeq1.length() * 100;

        // Print results
        System.out.println("Alignment Result (Affine Gap Penalty):");
        System.out.println("Aligned Sequence 1: " + alignedSeq1);
        System.out.println("Aligned Sequence 2: " + alignedSeq2);
        System.out.println("Length: " + alignedSeq1.length());
        System.out.printf("Identity: %d/%d (%.2f%%)%n", identities, alignedSeq1.length(), similarity);
        System.out.printf("Gaps: %d/%d (%.2f%%)%n", gaps, alignedSeq1.length(), gapPercent);
        System.out.println("Score: " + M[n][m]);
    }
    */
    //This function takes into consideration the gap scores + the new matrix EDNAFULL
    public void align_using_affine_gap_penalties() {
        int n = sequence1.length();
        int m = sequence2.length();

        // Initialize matrices
        double[][] M = new double[n + 1][m + 1]; // Match/mismatch matrix
        double[][] X = new double[n + 1][m + 1]; // Gap in sequence2
        double[][] Y = new double[n + 1][m + 1]; // Gap in sequence1
        int[][] tracebackMatrix = new int[n + 1][m + 1]; // To trace the alignment path

        // Initialize the matrices
        M[0][0] = 0;
        for (int i = 1; i <= n; i++) {
            M[i][0] = Double.NEGATIVE_INFINITY;
            X[i][0] = gap_open_penalty + (i - 1) * gap_extension_penalty;
            Y[i][0] = Double.NEGATIVE_INFINITY;
            tracebackMatrix[i][0] = 1; // Gap in sequence2
        }
        for (int j = 1; j <= m; j++) {
            M[0][j] = Double.NEGATIVE_INFINITY;
            X[0][j] = Double.NEGATIVE_INFINITY;
            Y[0][j] = gap_open_penalty + (j - 1) * gap_extension_penalty;
            tracebackMatrix[0][j] = 2; // Gap in sequence1
        }

        // Fill matrices
        for (int i = 1; i <= n; i++) {
            char char1 = sequence1.charAt(i - 1);
            for (int j = 1; j <= m; j++) {
                char char2 = sequence2.charAt(j - 1);

                // Calculate scores for M, X, and Y
                double scoreDiagonal = M[i - 1][j - 1] + getSubstitutionScore(char1, char2);
                double gapInSeq2 = Math.max(M[i - 1][j] + gap_open_penalty, X[i - 1][j] + gap_extension_penalty);
                double gapInSeq1 = Math.max(M[i][j - 1] + gap_open_penalty, Y[i][j - 1] + gap_extension_penalty);

                M[i][j] = Math.max(scoreDiagonal, Math.max(gapInSeq2, gapInSeq1));
                X[i][j] = gapInSeq2;
                Y[i][j] = gapInSeq1;

                // Update traceback matrix
                if (M[i][j] == scoreDiagonal) {
                    tracebackMatrix[i][j] = 0; // Match/Mismatch
                } else if (M[i][j] == gapInSeq2) {
                    tracebackMatrix[i][j] = 1; // Gap in sequence2
                } else {
                    tracebackMatrix[i][j] = 2; // Gap in sequence1
                }
            }
        }

        // Traceback to build aligned sequences
        StringBuilder alignedSeq1 = new StringBuilder();
        StringBuilder alignedSeq2 = new StringBuilder();
        int i = n, j = m;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && tracebackMatrix[i][j] == 0) {
                alignedSeq1.append(sequence1.charAt(i - 1));
                alignedSeq2.append(sequence2.charAt(j - 1));
                i--;
                j--;
            } else if (i > 0 && tracebackMatrix[i][j] == 1) {
                alignedSeq1.append(sequence1.charAt(i - 1));
                alignedSeq2.append('-');
                i--;
            } else {
                alignedSeq1.append('-');
                alignedSeq2.append(sequence2.charAt(j - 1));
                j--;
            }
        }

        alignedSeq1.reverse();
        alignedSeq2.reverse();

        // Calculate alignment statistics
        int identities = 0, gaps = 0;
        for (int k = 0; k < alignedSeq1.length(); k++) {
            if (alignedSeq1.charAt(k) == alignedSeq2.charAt(k)) {
                identities++;
            }
            if (alignedSeq1.charAt(k) == '-' || alignedSeq2.charAt(k) == '-') {
                gaps++;
            }
        }

        double similarity = (double) identities / alignedSeq1.length() * 100;
        double gapPercent = (double) gaps / alignedSeq1.length() * 100;

        // Print results
        System.out.println("Alignment Result (Affine Gap Penalty):");
        System.out.println("Aligned Sequence 1: " + alignedSeq1);
        System.out.println("Aligned Sequence 2: " + alignedSeq2);
        System.out.println("Length: " + alignedSeq1.length());
        System.out.printf("Identity: %d/%d (%.2f%%)%n", identities, alignedSeq1.length(), similarity);
        System.out.printf("Gaps: %d/%d (%.2f%%)%n", gaps, alignedSeq1.length(), gapPercent);
        System.out.println("Score: " + M[n][m]);
    }






    /////helper functions////////////////////////////////////////////////////////////////////////////////////////
    private double getOrDefault(Scanner scanner, double defaultValue) {
        String input = scanner.nextLine().trim();
        return input.isEmpty() ? defaultValue : Double.parseDouble(input);
    }

    public static String readFileContent(File file) throws IOException {
        return new String(Files.readAllBytes(file.toPath())).replaceAll("\\s", "");
    }

    public void info() {
        System.out.println("\nGalaxy tool needle: Needleman-Wunsch for Global Sequence Alignment");
        System.out.println("Sequence 1: " + sequence1);
        System.out.println("Sequence 2: " + sequence2);
        System.out.println("Insertion score = " + insertion);
        System.out.println("Deletion score = " + deletion);
        System.out.println("Mismatch score = " + mismatch);
        System.out.println("Match score = " + match);
        System.out.println("Gap open penalty = " + gap_open_penalty);
        System.out.println("Gap extension penalty = " + gap_extension_penalty);
    }





/////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public static void main(String[] args) {
        try {
            JFileChooser fileChooser = new JFileChooser();

            System.out.println("Choose the .txt file containing sequence 1");
            if (fileChooser.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) {
                System.out.println("File selection canceled.");
                return;
            }
            String sequence1 = readFileContent(fileChooser.getSelectedFile());

            System.out.println("Choose the .txt file containing sequence 2");
            if (fileChooser.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) {
                System.out.println("File selection canceled.");
                return;
            }
            String sequence2 = readFileContent(fileChooser.getSelectedFile());

            NeedlemanWunsch nw = new NeedlemanWunsch(sequence1, sequence2);
            nw.info();
            Scanner sc = new Scanner(System.in);
            System.out.println("There are two implemented versions for needleman-wunsch:\n1- Regular alignment (no affine gap penalty)\n2- Alignment with Affine Gap penalty\nWhich one would you rather run (Enter either 1 or 2)?");
            int choice = sc.nextInt();
            if (choice == 1) {
                nw.align();
            } else if (choice== 2){
                nw.align_using_affine_gap_penalties();
            } else {
                System.out.println("Wrong input");
            }
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }
    }

}



