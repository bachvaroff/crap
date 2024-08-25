	TITLE	$fasm

	.386
DGROUP	GROUP	CONST, _BSS, _DATA
EXTRN	__fltused:NEAR
PUBLIC  _fsincos
_DATA	SEGMENT  DWORD USE32 PUBLIC 'DATA'
_DATA      ENDS
_BSS	SEGMENT  DWORD USE32 PUBLIC 'BSS'
_BSS      ENDS
CONST	SEGMENT  DWORD USE32 PUBLIC 'CONST'
CONST      ENDS
_TEXT	SEGMENT  DWORD USE32 PUBLIC 'CODE'
	ASSUME   CS: _TEXT, DS: DGROUP, SS: DGROUP, ES: DGROUP

_fsincos	PROC NEAR
	fld	QWORD PTR [esp+4]
	fsincos
	mov	eax, DWORD PTR [esp+12]
	fstp	QWORD PTR [eax]
	mov	eax, DWORD PTR [esp+16]
	fstp	QWORD PTR [eax]
	ret
_fsincos	ENDP

_TEXT	ENDS
END
