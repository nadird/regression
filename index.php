
<?php
$text = "<h2>Notes :</h2>
<ul>
    <li>The entropy score for an alignment (a set of aligned sequences) is a measure of the variability at each position in the alignment. It quantifies the uncertainty or disorder in the sequence data. In the context of multiple sequence alignments (MSA), entropy is often used to evaluate conservation levels across aligned positions. A low entropy value indicates a highly conserved position where most sequences share the same character, suggesting functional or structural importance. Conversely, a high entropy value signifies greater variability, meaning the position is less conserved. This metric helps assess the reliability and significance of an alignment, particularly in evolutionary and structural biology studies. The function used is the one from the pymsa library in Python.</li>
    <li>In progress bars, we consider entropy values between 0 and 7000, which cover most cases (until 102 sequences by file and sequence lengths until 1200).</li>
    <li>Some examples of unaligned sequences and their reference aligned sequences from the HOMSTARD benchmark (233 pairs of files) are available at this link:  <a href='https://zenodo.org/records/16759607/files/HOMSTARD.zip?download=1'>Click here</a></li>
	

	
	
</ul>

<h2>FASTA File Criteria:</h2>
<h3>Unaligned Sequences</h3>
<ul>
    <li>Valid FASTA format: Each sequence starts with > followed by an identifier.</li>
    <li>No gaps (- or .).</li>
    <li>Only standard amino acids.</li>
    <li>Number of sequences: 3 - 102.</li>
    <li>Maximum sequence length: 1200.</li>
</ul>

<h3>Aligned Sequences</h3>
<ul>
    <li>Valid FASTA format: Each sequence starts with > followed by an identifier.</li>
    <li>Aligned sequences: All sequences must have the same length.</li>
    <li>Only standard amino acids and gaps (-).</li>
    <li>Number of sequences: 3 - 102.</li>
    <li>All sequences must have the same length.</li>
</ul>";
?>


<?php
session_start(); // Démarrer la session

// Fonction de validation basique du format FASTA
function isValidFasta1($content) {
    // Define the 20 standard amino acids
    $validAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    // Split content into lines
    $lines = explode("\n", trim($content));

    $sequences = [];
    $currentSequence = "";
    $sequenceCount = 0;

    foreach ($lines as $line) {
        $line = trim($line);
        
        if (empty($line)) continue; // Ignore empty lines
        
        if ($line[0] === '>') { // Detect a new sequence
            if (!empty($currentSequence)) {
                $sequences[] = $currentSequence;
                $currentSequence = "";
            }
            $sequenceCount++;
        } else {
            // Check that the sequence contains only valid amino acids
            if (strspn($line, $validAminoAcids) !== strlen($line)) {
                return "Error in unaligned sequences : Sequence contains invalid characters.";
            }
            $currentSequence .= $line;
        }
    }

    // Add the last sequence if it exists
    if (!empty($currentSequence)) {
        $sequences[] = $currentSequence;
    }

    // Check the number of sequences
    if ($sequenceCount < 3 || $sequenceCount > 102) {
        return "Error in unaligned sequences : The file must contain between 3 and 102 sequences (It contains $sequenceCount).";
    }

    // Check the length of each sequence
    foreach ($sequences as $seq) {
        if (strlen($seq) > 1200) {
            return "Error in unaligned sequences : A sequence exceeds the 1200-character limit.";
        }
    }

    return true; // Valid FASTA file
}



function isValidFasta2($content) {
    // Define allowed characters: 20 standard amino acids + gap ('-')
    $validCharacters = "ACDEFGHIKLMNPQRSTVWY-";

    // Split content into lines
    $lines = explode("\n", trim($content));

    $sequences = [];
    $currentSequence = "";
    $sequenceCount = 0;
    $sequenceLength = null;

    foreach ($lines as $line) {
        $line = trim($line);

        if (empty($line)) continue; // Ignore empty lines

        if ($line[0] === '>') { // New sequence identifier
            if (!empty($currentSequence)) {
                $sequences[] = $currentSequence;
                $currentSequence = "";
            }
            $sequenceCount++;
        } else {
			$line = strtoupper($line);   // Convert sequence to uppercase
            // Check if the sequence contains only valid characters
            if (strspn($line, $validCharacters) !== strlen($line)) {
                return "Error in aligned sequences : Sequence contains invalid characters.";
            }
            $currentSequence .= $line;
        }
    }

    // Add the last sequence if it exists
    if (!empty($currentSequence)) {
        $sequences[] = $currentSequence;
    }

    // Check the number of sequences
    if ($sequenceCount < 3 || $sequenceCount > 102) {
        return "Error in aligned sequences : The file must contain between 3 and 102 sequences.";
    }

    // Check if all sequences have the same length
    $sequenceLength = strlen($sequences[0]); // Reference length
    foreach ($sequences as $seq) {
        if (strlen($seq) !== $sequenceLength) {
            return "Error in aligned sequences : All sequences must have the same length.";
        }
    }

    return true; // Valid aligned FASTA file
}



// Définition des fonctions
function process($unalignedsequences) {
	$sequences = $unalignedsequences;
	$sequences2 = preg_split("/\R+/", trim($sequences)); // Utilise \R pour tout type de saut de ligne
	$sequences2 = array_values(array_filter($sequences2, 'strlen')); // Supprime les lignes vides et réindexe
	
	$file = "sequencesUnaln.tmp";
	$result = file_put_contents($file, implode("\n", $sequences2));
	
	$output = shell_exec("/home/nour_belghobsi/ml_env/bin/python predict.py " . escapeshellarg($file) . " 2>&1");
	//$_SESSION["error"] = $output;
	$sub_str = "CodeBegin";

	$output = strstr($output, $sub_str, false);
	$output = str_replace($sub_str, "", $output);

	// Split the output by comma to get individual results
		$results = explode(",", $output);
		$results[0]=trim($results[0]," []\n\r\t\v\x00");
		// Encodage des résultats pour éviter les erreurs dans l'URL
		//urlencode($results[0]);
		return $results[0];
		//ob_clean(); // Efface le buffer de sortie avant d'envoyer les headers
		//header("Location: index.php?resultA=$resultA");

}

function calculate($alignedsequences) {
	$sequences = $alignedsequences;
	$sequences2 = preg_split("/\R+/", trim($sequences)); // Utilise \R pour tout type de saut de ligne
	$sequences2 = array_values(array_filter($sequences2, 'strlen')); // Supprime les lignes vides et réindexe
	

	$file = "sequencesAlign.tmp";
	file_put_contents($file, implode("\n", $sequences2)); 
	$output = shell_exec("/home/nour_belghobsi/ml_env/bin/python msascorecalculate.py " . escapeshellarg($file));
	// Split the output by comma to get individual results
	$results = explode(",", $output);
	$results[0]=trim($results[0]," []\n\r\t\v\x00");
    //header("Location: index.php?score1=$results[0]");
	return $results[0];
}


// Vérifier quel bouton a été cliqué
if ($_SERVER["REQUEST_METHOD"] == "POST") {
	
	
    if (isset($_POST['process']) ) {
        $content = trim($_POST["fastaContent1"]);
		
		$validationResult = isValidFasta1($content);
		if ($validationResult === true) {
            $_SESSION["fastaContent1"] = htmlspecialchars($content);
            $_SESSION["result_A"] = process( htmlspecialchars_decode($_SESSION["fastaContent1"]));
		} else {
			$_SESSION["fastaContent1"] = htmlspecialchars($content);
			$_SESSION["result_A"]  ="";
			 $_SESSION["error"] = $validationResult; // Display the error message
		}
		
        header("Location: " . $_SERVER['PHP_SELF']);
        exit();
    } elseif (isset($_POST['calculate']) ) {
        $content = trim($_POST["fastaContent2"]);
		
		$validationResult = isValidFasta2($content);
		if ($validationResult === true) {
            $_SESSION["fastaContent2"] = htmlspecialchars($content);
            $_SESSION["score_1"] = calculate( htmlspecialchars_decode($_SESSION["fastaContent2"]));
		} else {
			$_SESSION["fastaContent2"] = htmlspecialchars($content);
			$_SESSION["score_1"] ="";
			 $_SESSION["error"] = $validationResult; // Display the error message
		}
		
        header("Location: " . $_SERVER['PHP_SELF']);
        exit();
    }
}



// Récupération des valeurs en session
$result_A = isset($_SESSION["result_A"]) ? $_SESSION["result_A"] : "";
$score_1 = isset($_SESSION["score_1"]) ? $_SESSION["score_1"] : "";
$error = isset($_SESSION["error"]) ? $_SESSION["error"] : "";
unset($_SESSION["error"]); // Supprimer le message d'erreur après affichage
?>

<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Entropy score prediction</title>
	
	 <style>
        body {
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f9f9f9;
            color: #333;
        }
        .container {
            max-width: 800px;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
        }
        h2 {
            color: #0073e6;
        }
        h3 {
            color: #005bb5;
        }
        ul {
            line-height: 1.6;
        }
        a {
            color: #0073e6;
            text-decoration: none;
            font-weight: bold;
        }
        a:hover {
            text-decoration: underline;
        }
		progress {
        width: 100%; /* Full width of the container */
        height: 20px; /* Adjust height if needed */
		}		
    </style>
	
    <script>
		function readFile1(input, targetTextbox) {
			if (input.files && input.files[0]) {
				var reader = new FileReader();
				reader.onload = function (e) {
					var content = e.target.result; // Lire le fichier tel quel
					document.getElementById(targetTextbox).value = content;
				};
				reader.readAsText(input.files[0]);
				document.getElementById("result_a").value = "";
				document.getElementById("myErrorLabel").textContent = "";
			}
		}
		
		function readFile2(input, targetTextbox) {
			if (input.files && input.files[0]) {
				var reader = new FileReader();
				reader.onload = function (e) {
					var content = e.target.result; // Lire le fichier tel quel
					document.getElementById(targetTextbox).value = content;
				};
				reader.readAsText(input.files[0]);
				document.getElementById("scr_1").value = "";
				document.getElementById("myErrorLabel").textContent = "";				
			}
		}
    </script>
</head>
<body>

<style>
    input[type="file"] {
        color: transparent;
    }
</style>

    <p style="color: red;"><?php echo "<label id='myErrorLabel' for='myInput'>$error</label>"; ?></p>
		
    <h2>Predicting the quality (Enropy score) of the best expected alignment from unaligned protein sequences :</h2>
	<label for="aa">- (The prediction is generated using artificial intelligence techniques, leveraging Artifficial Neural Networks (ANNs) to process the input data.)</label> 
	<br><br>
	

    <form method="post">
        <label for="fasta1">Upload unaligned sequences (Fasta file) :</label><br>
        <input type="file" id="fasta1" name="fasta1" onchange="readFile1(this, 'fastaContent1')"><br>
        <textarea id="fastaContent1" name="fastaContent1" rows="10" cols="100" readonly><?php echo isset($_SESSION["fastaContent1"]) ? $_SESSION["fastaContent1"] : ''; ?></textarea><br>
        <button type="submit" name="process">Predict</button>
    </form>

    <label for="result_a">Predicted entropy score:</label>
    <input type="text" id="result_a" value="<?php echo htmlspecialchars($result_A); ?>" disabled>
	<br>
	<label for="grandeur">------- :</label>
    <progress id="grandeur" value="<?php echo htmlspecialchars($result_A); ?>" max="6000"></progress>
	
	<br><br><br>
	
    <h2>Calulate the quality (Enropy score) of a given alignment :</h2>
    <form method="post">
        <label for="fasta2">Upload aligned sequences (Fasta file) :</label><br>
        <input type="file" id="fasta2" name="fasta2" onchange="readFile2(this, 'fastaContent2')"><br>
        <textarea id="fastaContent2" name="fastaContent2" rows="10" cols="100"  readonly><?php echo isset($_SESSION["fastaContent2"]) ? $_SESSION["fastaContent2"] : ''; ?></textarea><br>
        <button type="submit" name="calculate">Calculate</button>
    </form>

    <label for="scr_1">Calculated entropy score :</label>
    <input type="text" id="scr_1" value="<?php echo htmlspecialchars($score_1); ?>" disabled>
	<br>
	<label for="grandeur">-------  :</label>
    <progress id="grandeur" value="<?php echo htmlspecialchars($score_1); ?>" max="6000"></progress>
	
	<br><br>
    <div class="container">
        <?php echo $text; ?>
    </div>
</body>
<footer>
    <p>&copy; <?php echo "2025"; ?> - Mohamed Skander Daas</p>
    <p>Last modified: <?php echo "August 7, 2025"; ?></p>
</footer>
</html>
