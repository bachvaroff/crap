TITLE $fasm

.386
.387

DGROUP GROUP CONST, _BSS, _DATA

EXTRN __fltused:NEAR

_DATA SEGMENT DWORD USE32 PUBLIC 'DATA'
_DATA ENDS

_BSS SEGMENT DWORD USE32 PUBLIC 'BSS'
_BSS ENDS

CONST SEGMENT DWORD USE32 PUBLIC 'CONST'
CONST ENDS

_TEXT SEGMENT DWORD USE32 PUBLIC 'CODE'
ASSUME CS: _TEXT, DS: DGROUP, SS: DGROUP, ES: DGROUP

_fsincos PROC NEAR
	fld qword ptr [esp + 4]
	fsincos
	mov eax, dword ptr [esp + 12]
	fstp qword ptr [eax]
	mov eax, dword ptr [esp + 16]
	fstp qword ptr [eax]
	ret
_fsincos ENDP

_TEXT ENDS

PUBLIC _fsincos

END

