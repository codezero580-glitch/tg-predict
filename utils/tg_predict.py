import numpy as np

from pathlib import Path
from tensorflow.keras.models import load_model

BASE_DIR = Path(__file__).resolve().parent.parent

model_path = BASE_DIR / "model" / "best_model_30_25.keras"

best_model = load_model(model_path)


def TokenizeSMILESasInt(smiles, token2idx):
    char_list = list(smiles)
    tokenized_arr = []

    for char in char_list:
        if char in token2idx:
            tokenized_arr.append(token2idx[char])
        else:
            print(f"Karakter tidak dikenal: {repr(char)}")

    return np.array(tokenized_arr, dtype=np.int32)


def TokenPolymerInt(SMI_list, token2idx, maxLen):
    Nsamples = len(SMI_list)
    Tokened_polymer = np.zeros(
        (Nsamples, maxLen),
        dtype=np.int32
    )

    for i in range(Nsamples):
        ismi = SMI_list[i]
        iarr = TokenizeSMILESasInt(ismi, token2idx)

        if iarr.shape[0] < maxLen:
            Tokened_polymer[i, :iarr.shape[0]] = iarr
        else:
            Tokened_polymer[i, :] = iarr[:maxLen]

    return Tokened_polymer


def predict_tg(psmiles, token2idx, maxLen):
    X = TokenPolymerInt(
        [psmiles],
        token2idx,
        maxLen
    )

    # Prediksi satu sampel
    y_pred = best_model.predict(X, verbose=0)

    return float(np.asarray(y_pred).reshape(-1)[0])
